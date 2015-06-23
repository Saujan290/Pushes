use Modern::Perl;
use HackaMol;
use Data::Dumper;

my $mol = HackaMol->new->read_file_mol(shift);

my @atoms1 = map{$mol->get_atoms($_)} (4,5,11,10);
my ($dihe1) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (4,5,11,10));
my ($dihe2) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (4,5,11,10));

$dihe1->dihe_fc(3.5);
$dihe1->dihe_mult(2);
$dihe2->dihe_fc(0.6);
$dihe2->dihe_mult(3);

my @atoms2 = map{$mol->get_atoms($_)} (2,1,4,5);
my ($dihe3) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (2,1,4,5));

$dihe3->dihe_fc(2);
$dihe3->dihe_mult(3);

my @atoms3 = map{$mol->get_atoms($_)} (11,10,7,8);
my ($dihe4) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (11,10,7,8));

$dihe4->dihe_fc(2);
$dihe4->dihe_mult(3);

my @atoms4 = map{$mol->get_atoms($_)} (1,4,5,11);
my ($dihe5) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (1,4,5,11));

$dihe5->dihe_fc(1);
$dihe5->dihe_mult(3);

my @atoms5 = map{$mol->get_atoms($_)} (5,11,10,7);
my ($dihe6) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (5,11,10,7));

$dihe6->dihe_fc(1);
$dihe6->dihe_mult(3);


my @biggens;
my @arry;
foreach my $t (0 .. $mol->tmax){
  $mol->t($t);
  my $ener = $dihe1->torsion_energy + $dihe2->torsion_energy + $dihe3->torsion_energy + $dihe4->torsion_energy + $dihe5->torsion_energy + $dihe6->torsion_energy;
  #my $dse = 2.0 * (1+ cos(3 * $dihe->dihe_rad));
  #push @biggens, $ener if ($ener >2 and $ener <3);
  #  $dihe->dihe_deg, "\t", $ener ; 
  push @arry, sprintf("%10.2f %10.2f\n", $dihe1->dihe_deg,$ener); 
};

print foreach @arry;



#my $sum += $sum + $_ ;
#};
#say $sum;
#say 'shit ' , $_ foreach @biggens;
