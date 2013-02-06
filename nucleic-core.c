/*

  File: "nucleic.c"

  Authors: Marcel Turcotte and Marc Feeley

  Last Modified On: Sat Mar 12 13:05:05 1994

  This program is described in the paper:

    "Using Multilisp for Solving Constraint Satisfaction Problems: an
    Application to Nucleic Acid 3D Structure Determination" published in
    the journal "Lisp and Symbolic Computation".

  This code has been hacked into multiple files (by Paul Kelly) to 
  get it to build under simplescalar, which has a problem with 
  relocation overflow in very large assembly code files.

  Further details of the benchmark can be found in this article:

    "Benchmarking implementations of functional languages with ‘Pseudoknot’, a     float-intensive benchmark
    Pieter H. Hartel, Marc Feeley, Martin Alt, Lennart Augustsson, 
    Peter Baumann, Marcel Beemster, Emmanuel Chailloux, Christine H. Flood, 
    Wolfgang Grieskamp, John H. G. Van Groningen, Kevin Hammond, Bogumil 
    Hausman, Melody Y. Ivory, Richard E. Jones, Jasper Kamperman, Peter Lee, 
    Xavier Leroy, Rafael D. Lins, Sandra Loosemore, Niklas Röjemo, Manuel 
    Serrano, Jean-Pierre Talpin, Jon Thackray, Stephen Thomas, Pum Walters, 
    Pierre Weis and Peter Wentworth.
    Journal of Functional Programming / Volume 6 / Issue 04 / July 1996, 
    pp 621-655"
    http://dx.doi.org/10.1017/S0956796800001891
    Available from http://doc.utwente.nl/55704/1/JFP_pseudoknotI.pdf
*/

#include <stdio.h>
#include <math.h>
#include "/tmp/nucleic.h"

void pt_sub( p1, p2, p3 )
pt p1, p2, p3;
{
  p3->x = p1->x - p2->x;
  p3->y = p1->y - p2->y;
  p3->z = p1->z - p2->z;
}

FLOAT pt_dist( p1, p2 )
pt p1, p2;
{
  FLOAT dx = p1->x - p2->x, dy = p1->y - p2->y, dz = p1->z - p2->z;
  return sqrt( dx*dx + dy*dy + dz*dz );
}

FLOAT pt_phi( p )
pt p;
{
  FLOAT x = p->x, y = p->y, z = p->z;
  FLOAT b = atan2( x, z );
  return atan2( cos(b)*z + sin(b)*x, y );
}

FLOAT pt_theta( p )
pt p;
{
  return atan2( p->x, p->z );
}

/*---------------------------------------------------------------------------*/

/* COORDINATE TRANSFORMATIONS */

void tfo_apply( t, p1, p2 )
tfo t;
pt p1, p2;
{
  FLOAT x = p1->x, y = p1->y, z = p1->z;
  p2->x = x*t->m[0][0] + y*t->m[1][0] + z*t->m[2][0] + t->m[3][0];
  p2->y = x*t->m[0][1] + y*t->m[1][1] + z*t->m[2][1] + t->m[3][1];
  p2->z = x*t->m[0][2] + y*t->m[1][2] + z*t->m[2][2] + t->m[3][2];
}

void tfo_combine( a, b, t )
tfo a, b, t;
{
  FLOAT x, y, z;

  x = a->m[0][0];  y = a->m[0][1];  z = a->m[0][2];
  t->m[0][0] = x*b->m[0][0] + y*b->m[1][0] + z*b->m[2][0];
  t->m[0][1] = x*b->m[0][1] + y*b->m[1][1] + z*b->m[2][1];
  t->m[0][2] = x*b->m[0][2] + y*b->m[1][2] + z*b->m[2][2];

  x = a->m[1][0];  y = a->m[1][1];  z = a->m[1][2];
  t->m[1][0] = x*b->m[0][0] + y*b->m[1][0] + z*b->m[2][0];
  t->m[1][1] = x*b->m[0][1] + y*b->m[1][1] + z*b->m[2][1];
  t->m[1][2] = x*b->m[0][2] + y*b->m[1][2] + z*b->m[2][2];

  x = a->m[2][0];  y = a->m[2][1];  z = a->m[2][2];
  t->m[2][0] = x*b->m[0][0] + y*b->m[1][0] + z*b->m[2][0];
  t->m[2][1] = x*b->m[0][1] + y*b->m[1][1] + z*b->m[2][1];
  t->m[2][2] = x*b->m[0][2] + y*b->m[1][2] + z*b->m[2][2];

  x = a->m[3][0];  y = a->m[3][1];  z = a->m[3][2];
  t->m[3][0] = x*b->m[0][0] + y*b->m[1][0] + z*b->m[2][0] + b->m[3][0];
  t->m[3][1] = x*b->m[0][1] + y*b->m[1][1] + z*b->m[2][1] + b->m[3][1];
  t->m[3][2] = x*b->m[0][2] + y*b->m[1][2] + z*b->m[2][2] + b->m[3][2];
}

void tfo_inv_ortho( t1, t2 )
tfo t1, t2;
{
  FLOAT tx = t1->m[3][0], ty = t1->m[3][1], tz = t1->m[3][2];
  t2->m[0][0] = t1->m[0][0];
  t2->m[0][1] = t1->m[1][0];
  t2->m[0][2] = t1->m[2][0];
  t2->m[1][0] = t1->m[0][1];
  t2->m[1][1] = t1->m[1][1];
  t2->m[1][2] = t1->m[2][1];
  t2->m[2][0] = t1->m[0][2];
  t2->m[2][1] = t1->m[1][2];
  t2->m[2][2] = t1->m[2][2];
  t2->m[3][0] = -(tx*t1->m[0][0] + ty*t1->m[0][1] + tz*t1->m[0][2]);
  t2->m[3][1] = -(tx*t1->m[1][0] + ty*t1->m[1][1] + tz*t1->m[1][2]);
  t2->m[3][2] = -(tx*t1->m[2][0] + ty*t1->m[2][1] + tz*t1->m[2][2]);
}

void tfo_align( p1, p2, p3, t )
pt p1, p2, p3;
tfo t;
{
  FLOAT x1 = p1->x, y1 = p1->y, z1 = p1->z;
  FLOAT x3 = p3->x, y3 = p3->y, z3 = p3->z;
  FLOAT x31 = x3-x1, y31 = y3-y1, z31 = z3-z1;
  struct struct_pt rotpy;
  pt_sub( p2, p1, &rotpy );
  {
    FLOAT phi = pt_phi( &rotpy ), theta = pt_theta( &rotpy );
    FLOAT sinp = sin(phi), sint = sin(theta);
    FLOAT cosp = cos(phi), cost = cos(theta);
    FLOAT sinpsint = sinp*sint;
    FLOAT sinpcost = sinp*cost;
    FLOAT cospsint = cosp*sint;
    FLOAT cospcost = cosp*cost;
    struct struct_pt rotpz;
    make_pt( cost*x31 - sint*z31,
             sinpsint*x31 + cosp*y31 + sinpcost*z31,
             cospsint*x31 - sinp*y31 + cospcost*z31,
             &rotpz );
    {
      FLOAT rho = pt_theta( &rotpz );
      FLOAT cosr = cos( rho ), sinr = sin( rho );
      FLOAT x = -x1*cost + z1*sint;
      FLOAT y = -x1*sinpsint - y1*cosp - z1*sinpcost;
      FLOAT z = -x1*cospsint + y1*sinp - z1*cospcost;
      t->m[0][0] = cost*cosr - cospsint*sinr;
      t->m[0][1] = sinpsint;
      t->m[0][2] = cost*sinr + cospsint*cosr;
      t->m[1][0] = sinp*sinr;
      t->m[1][1] = cosp;
      t->m[1][2] = -sinp*cosr;
      t->m[2][0] = -sint*cosr - cospcost*sinr;
      t->m[2][1] = sinpcost;
      t->m[2][2] = -sint*sinr + cospcost*cosr;
      t->m[3][0] = x*cosr - z*sinr;
      t->m[3][1] = y;
      t->m[3][2] = x*sinr + z*cosr;
    }
  }
}



/*---------------------------------------------------------------------------*/

/* SEARCH */

void (*seq)();
int (*constraint)();
sol_list solutions, last_solution;

void search( problem_domains, problem_constraint )
void (*problem_domains)();
int (*problem_constraint)();
{
  seq = problem_domains;
  constraint = problem_constraint;
  solutions = NULL;
  last_solution = NULL;
  partial_inst_length = 0;
  (*seq)();
}

void try( i, t, n )
int i;
tfo t;
nuc n;
{
  var v = &partial_inst[partial_inst_length];
  v->i = i;
  v->t = t;
  v->n = n;
  if ((*constraint)( v ))
  {
    partial_inst_length++;
/*    printf("partial_inst_length=%d\n", partial_inst_length); */
    (*seq)();
    partial_inst_length--;
  }
}

void found_solution()
{
  int i;
  sol_list sols;
  var_list sol;
  sol = (var_list)malloc( (partial_inst_length+1)*sizeof(struct struct_var) );
  if (sol==NULL) { printf( "memory overflow\n" ); exit(1); }
  for (i=0; i<partial_inst_length; i++)
  {
    tfo t;
    t = (tfo)malloc( sizeof(struct struct_tfo) );
    if (t==NULL) { printf( "memory overflow\n" ); exit(1); }
    *t = *partial_inst[i].t;
    sol[i].i = partial_inst[i].i;
    sol[i].t = t;
    sol[i].n = partial_inst[i].n;
  }
  sol[i].i = -1;
  sols = (sol_list)malloc( sizeof(struct struct_sol_list) );
  if (sols==NULL) { printf( "memory overflow\n" ); exit(1); }
  sols->sol = sol;
  sols->next = NULL;
  if (last_solution == NULL)
    solutions = sols;
  else
    last_solution->next = sols;
  last_solution = sols;
}


/* IO */

int atom_num;

void dump_pt( out, t, nuc_num, nuc_name, atom_name, p )
FILE *out;
tfo t;
int nuc_num;
char nuc_name, *atom_name;
pt p;
{
  struct struct_pt atom;
  tfo_apply( t, p, &atom );
  fprintf( out,
           "ATOM  %5d  %-4s %c   %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           atom_num, atom_name, nuc_name, nuc_num,
           atom.x, atom.y, atom.z, 1.0, 0.0 );
  atom_num++;
}

void dump_var( out, v )
FILE *out;
var v;
{
  int i = v->i;
  tfo t = v->t;
  nuc n = v->n;
  dump_pt( out, t, i, n->type, "P",   &n->p    );
  dump_pt( out, t, i, n->type, "O1P", &n->o1p  );
  dump_pt( out, t, i, n->type, "O2P", &n->o2p  );
  dump_pt( out, t, i, n->type, "O5*", &n->o5_  );
  dump_pt( out, t, i, n->type, "C5*", &n->c5_  );
  dump_pt( out, t, i, n->type, "H5*", &n->h5_  );
  dump_pt( out, t, i, n->type, "H5**",&n->h5__ );
  dump_pt( out, t, i, n->type, "C4*", &n->c4_  );
  dump_pt( out, t, i, n->type, "H4*", &n->h4_  );
  dump_pt( out, t, i, n->type, "O4*", &n->o4_  );
  dump_pt( out, t, i, n->type, "C1*", &n->c1_  );
  dump_pt( out, t, i, n->type, "H1*", &n->h1_  );
  dump_pt( out, t, i, n->type, "C2*", &n->c2_  );
  dump_pt( out, t, i, n->type, "H2**",&n->h2__ );
  dump_pt( out, t, i, n->type, "O2*", &n->o2_  );
  dump_pt( out, t, i, n->type, "H2*", &n->h2_  );
  dump_pt( out, t, i, n->type, "C3*", &n->c3_  );
  dump_pt( out, t, i, n->type, "H3*", &n->h3_  );
  dump_pt( out, t, i, n->type, "O3*", &n->o3_  );
  dump_pt( out, t, i, n->type, "N1",  &n->n1   );
  dump_pt( out, t, i, n->type, "N3",  &n->n3   );
  dump_pt( out, t, i, n->type, "C2",  &n->c2   );
  dump_pt( out, t, i, n->type, "C4",  &n->c4   );
  dump_pt( out, t, i, n->type, "C5",  &n->c5   );
  dump_pt( out, t, i, n->type, "C6",  &n->c6   );

  switch (n->type)
  {
    case 'A':
      dump_pt( out, t, i, n->type, "N6",  &n->_.A.n6   );
      dump_pt( out, t, i, n->type, "N7",  &n->_.A.n7   );
      dump_pt( out, t, i, n->type, "N9",  &n->_.A.n9   );
      dump_pt( out, t, i, n->type, "C8",  &n->_.A.c8   );
      dump_pt( out, t, i, n->type, "H2",  &n->_.A.h2   );
      dump_pt( out, t, i, n->type, "H61", &n->_.A.h61  );
      dump_pt( out, t, i, n->type, "H62", &n->_.A.h62  );
      dump_pt( out, t, i, n->type, "H8",  &n->_.A.h8   );
      break;

    case 'C':
      dump_pt( out, t, i, n->type, "N4",  &n->_.C.n4   );
      dump_pt( out, t, i, n->type, "O2",  &n->_.C.o2   );
      dump_pt( out, t, i, n->type, "H41", &n->_.C.h41  );
      dump_pt( out, t, i, n->type, "H42", &n->_.C.h42  );
      dump_pt( out, t, i, n->type, "H5",  &n->_.C.h5   );
      dump_pt( out, t, i, n->type, "H6",  &n->_.C.h6   );
      break;

    case 'G':
      dump_pt( out, t, i, n->type, "N2",  &n->_.G.n2   );
      dump_pt( out, t, i, n->type, "N7",  &n->_.G.n7   );
      dump_pt( out, t, i, n->type, "N9",  &n->_.G.n9   );
      dump_pt( out, t, i, n->type, "C8",  &n->_.G.c8   );
      dump_pt( out, t, i, n->type, "O6",  &n->_.G.o6   );
      dump_pt( out, t, i, n->type, "H1",  &n->_.G.h1   );
      dump_pt( out, t, i, n->type, "H21", &n->_.G.h21  );
      dump_pt( out, t, i, n->type, "H22", &n->_.G.h22  );
      dump_pt( out, t, i, n->type, "H8",  &n->_.G.h8   );
      break;

    case 'U':
      dump_pt( out, t, i, n->type, "O2",  &n->_.U.o2   );
      dump_pt( out, t, i, n->type, "O4",  &n->_.U.o4   );
      dump_pt( out, t, i, n->type, "H3",  &n->_.U.h3   );
      dump_pt( out, t, i, n->type, "H5",  &n->_.U.h5   );
      dump_pt( out, t, i, n->type, "H6",  &n->_.U.h6   );
      break;
  }
}

void dump_solution( out, sol )
FILE *out;
var_list sol;
{
  atom_num = 1;
  while (sol->i >= 0)
  {
    dump_var( out, sol );
    sol++;
  }
}

int count_solutions( sl, save )
sol_list sl;
int save;
{
  int n;
  sol_list l;

  n = 0;
  l = sl;
  while (l!=NULL)
  {
    if (save)
    {
      char filename[30];
      FILE *out;
      sprintf( filename, "sol%04d.pdb", n+1 );
      out = fopen( filename, "w" );
      dump_solution( out, l->sol );
      fclose( out );
    }
    l = l->next;
    n++;
  }
  return n;
}

/*---------------------------------------------------------------------------*/

/* TESTING */

/*

  The problem to be solved is selected through the first command-line
  argument.  The second argument controls the dumping of the solutions
  to disk.  Assuming "nucleic" is the name of the binary,

    "nucleic 1" - anticodon     runs the anticodon problem
    "nucleic 4" - pseudoknot    runs the pseudoknot problem

  The result printed by the program should be 179 for anticodon and 50
  for pseudoknot.

  To save the solutions out to disk in PDB format use:

    nucleic 1 save
    nucleic 4 save

  The resulting files (named "solNNNN.pdb") can be viewed with the following
  public domain programs:

   XMOL -- Minnesota Supercomputer Center
   WHERE: ftp.msc.edu:pub/xmol   (no source)
   HOSTS: mips, rs6000, sgi, sun4

   RasMol -- University of Edinburgh
   WHERE: ftp.dcs.ed.ac.uk:pub/rasmol   (source)
   HOSTS: sun3, sun4, sun386i, hp9000, sequent, DEC alpha, IBM RS/6000
          and SGI, DEC and E&S, Linux, MS Windows

   tkrasmol (tcl/tk X toolkit rasmol) -- University of Berkeley
   WHERE: sprite.Berkeley.EDU:tkrasmol.tar.Z   (source)
   HOSTS: any UNIX-like system that approximates POSIX, BSD, or System V

   Kinemage -- University of California at Irvine
   WHERE: orion.oac.uci.edu:/protein/Kinemage   (no source)
   HOSTS: Windows, Macintosh

*/

/*int equal( str1, str2 )
char *str1, *str2;
{
  printf("checking command line argument %s == %s ?\n", str1, str2);
  while ((*str1!='\0') && (*str2!='\0') && (*str1==*str2)) { str1++; str2++; }
  printf("Done\n");
  return (*str1==*str2);
}*/

int main(argc, argv)
int argc;
char **argv;
{
  extern long clock();
  long lasttime;

  if (argc == 1)
  {
    printf( "Usage: %s 1|2|3|4 [save]\n", argv[0] );
    exit(1);
  }

  /* Argument parsing has been hacked horribly in an unsuccessful attempt to 
     overcome problems with the simplescalar compiler - Paul Kelly */

  switch(argv[1][0]) {
  case '1':
    {
      printf("anticodon\n");
      init_nucleotides();
      lasttime = clock();
      search( anticodon_domains, anticodon_constraint );
      printf("anticodon: %ld\n", clock()-lasttime);
      break;
    }
  case '2':
    {
      printf("anticodonV2\n");
      init_nucleotides();
      lasttime = clock();
      search( anticodonV2_domains, anticodon_constraint );
      printf("anticodonV2: %ld\n", clock()-lasttime);
      break;
    }
  case '3':
    {
      printf("anticodonV3");
      init_nucleotides();
      lasttime = clock();
      search( anticodonV3_domains, anticodon_constraint );
      printf("anticodonV3: %ld\n", clock()-lasttime);
      break;
    }
  case '4':
    {
      printf("pseudoknot\n");
      init_nucleotides();
      lasttime = clock();
      search( pseudoknot_domains, pseudoknot_constraint );
      printf("pseudoknot: %ld\n", clock()-lasttime);
      break;
    }
  default:
    {
      printf("argument not understood\n");
      exit(1);
    }
  }
  printf( "Number of solutions: %d\n", count_solutions( solutions, argc>2 ) );

  return 0;
}

/*---------------------------------------------------------------------------*/
