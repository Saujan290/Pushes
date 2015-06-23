use HackaMol;
use HackaMol::X::Orca;
use  HackaMol::X::NERF;
use Modern::Perl;
use YAML::XS qw(LoadFile DumpFile);
use POSIX;

# setup

my $t1 = time;

my $file = shift;
my $mol = HackaMol->new->read_file_mol($file);
say 'Done Reading';

my $root = $file;
$root =~ s/\.xyz//;
my $histo = LoadFile("xtals_best.txt");
my %histo = %{$histo};


#establish dihedrals
my $dihe1 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(1),
                                $mol->get_atoms(0),
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                             ],
);

my $dihe2 = HackaMol::Dihedral->new(
                               is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(0),
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                $mol->get_atoms(10),
                                             ],
);

my $dihe3 = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(3),
                                $mol->get_atoms(4),
                                $mol->get_atoms(10),
                                $mol->get_atoms(9),
                                             ],
);

my $dihe1p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                   atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(6),
                                $mol->get_atoms(9),
                                $mol->get_atoms(10),
                                             ],
);


my $dihe2p = HackaMol::Dihedral->new(
                                is_constrained => 1,
                                    atoms=> [
                                $mol->get_atoms(7),
                                $mol->get_atoms(9),
                                $mol->get_atoms(10),
                               $mol->get_atoms(4),
                                             ],
);



# push dihedrals
$mol->push_dihedrals($dihe1);
$mol->push_dihedrals($dihe2);
$mol->push_dihedrals($dihe3);
$mol->push_dihedrals($dihe1p);
$mol->push_dihedrals($dihe2p);

#checks
say $dihe1->dihe_deg;
say $dihe2->dihe_deg;
say $dihe1p->dihe_deg;
say $dihe2p->dihe_deg;
say $dihe3->dihe_deg;

#set up orca
my $orca2 = HackaMol::X::Orca->new(
      mol             => $mol,
      has_constraints => 1,      
      theory          => 'HF-3c',
      exe             => '/Users/chem_student/perl5/apps/orca_3_0_3_macosx_openmpi165/orca',
      scratch         => 'tmp',
);

# get ready for the foreach loop
$mol->t(0);
open(my $in, ">", "xtals_energies.txt") or die "couldn't open";
#my %fresults;
#my $partit = 5;
#my $partit2 = 1;

#my $nnbond = HackaMol::Bond->new(
#    atoms=> [
#    $mol->get_atoms(6),
#   $mol->get_atoms(0),
#  ]
#);


#say $nnbond->bond_length;
#my $nnlength = 0;

foreach my $dist (sort { $a <=> $b } keys %histo){
  foreach my $chi3 (sort{ $a <=> $b } keys %{$histo{$dist}}){
      my $t = $histo{$dist}{$chi3};
      my $qm_mol = cystine_protonate($mol, $t);
    #foreach my $t (@{$histo{$dist}{$chi3}}){
      #my @energies = (0);
      my @energies = $orca2->opt;
      my $mol2 = $orca2->load_trj;
      $mol2->print_xyz_ts([0 .. $mol2->tmax]);
      printf ("%10.3f %10.3f %10.2f\n", $dist, $chi3,   $mol2->get_energy($mol2->tmax)*627.51);
      printf $in ("%10.3f %10.3f %10.2f\n",$dist, $chi3,  $energies[-1]*627.51); 
      exit;
  } 
}

sub cystine_protonate {
  my $mol = shift;
  my $t   = shift;

  $mol->t($t);

  my $nerf = HackaMol::X::NERF->new;
  #build in two oxygen
  my @xyzs = map{$_->xyz} $mol->all_atoms;
  use Data::Dumper;
  say 'shit'; 
  say $_ foreach @xyzs;

  my $o1  = $nerf->extend_abc( $xyzs[3],  $xyzs[1], $xyzs[2],  1.3,  120, 180 );
  my $o2  = $nerf->extend_abc( $xyzs[9],  $xyzs[7], $xyzs[8],  1.3,  120, 180 );
  my $hmol = HackaMol::Molecule->new(atoms=>[$mol->all_atoms,
        map{ HackaMol::Atom->new(Z => 8, coords => [$_]) } ($o1,$o2)
      ]
  );
  
  @xyzs = map{$_->xyz} $hmol->all_atoms;
  my $h1  = $nerf->extend_abc( $xyzs[2],  $xyzs[1], $xyzs[0],  1.0,  120, 1 );
  my $h2  = $nerf->extend_abc( $xyzs[2],  $xyzs[1], $xyzs[0],  1.0,  120, 180 );
  my $h3  = $nerf->extend_abc( $xyzs[7],  $xyzs[5], $xyzs[6],  1.0,  120, 1 );
  my $h4  = $nerf->extend_abc( $xyzs[7],  $xyzs[5], $xyzs[6],  1.0,  120, 180 );
  my $h5  = $nerf->extend_abc( $xyzs[2],  $xyzs[0], $xyzs[1],  1.0,  109, 120);
  my $h6  = $nerf->extend_abc( $xyzs[8],  $xyzs[6], $xyzs[7],  1.0,  109, 90);
  my $h7  = $nerf->extend_abc( $xyzs[11],  $xyzs[5], $xyzs[4],  1.0,  109, -90);
  my $h8  = $nerf->extend_abc( $xyzs[11],  $xyzs[5], $xyzs[4],  1.0,  109, 90);
  my $h9  = $nerf->extend_abc( $xyzs[11], $xyzs[9], $xyzs[10], 1.0,  109, 270 );
  my $h10 = $nerf->extend_abc( $xyzs[11], $xyzs[9], $xyzs[10], 1.0,  109, 90);
  my $h11 = $nerf->extend_abc( $xyzs[3],  $xyzs[2], $xyzs[12], 1.0,  109, 180);
  my $h12 = $nerf->extend_abc( $xyzs[9],  $xyzs[8], $xyzs[13], 1.0,  109, 54);

  $hmol->push_atoms(
    map{ HackaMol::Atom->new(Z => 1, coords => [$_]) } ($h1,$h2,$h3,$h4,$h5,$h6,$h7,$h8,$h9,$h10,$h11,$h12)
  );
#  my $hmol = HackaMol::Molecule->new(atoms=>[$mol->all_atoms,
#      map{ HackaMol::Atom->new(Z => 8, coords => [$_]) } ($o1,$o2),
#      map{ HackaMol::Atom->new(Z => 1, coords => [$_]) } ($h1,$h2,$h3,$h4,$h5,$h6,$h7,$h8,$h9,$h10,$h11,$h12)
#      ]
#  ); 
 
  $hmol->print_xyz;
  exit;

  return ($hmol);
  
}
  
  

