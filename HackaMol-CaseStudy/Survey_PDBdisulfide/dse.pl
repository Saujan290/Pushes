use Modern::Perl;
use HackaMol;
use Data::Dumper;

my $mol = HackaMol->new->read_file_mol(shift);

my @atoms2 = map{$mol->get_atoms($_)} (1,6,7,8);
my ($dihe1) = HackaMol->new->build_dihedrals( map {$mol->get_atoms($_)} (1,6,7,8));

$dihe1->dihe_fc(2);
$dihe1->dihe_mult(3);

my @biggens;
my @arry;
foreach my $t (0 .. $mol->tmax){
  $mol->t($t);
  my $ener = $dihe1->torsion_energy;
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
