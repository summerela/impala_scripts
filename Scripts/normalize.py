__author__ = 'selasady'
'''Normalization algorithm based on University of Michigan algorithm based on perl script
 written for Kaviar by Terry Farrah'''







# sub trimVar {
#   # typically, the first element of @alleles is ref
#   my ($chrom, $frz, $start, @alleles) = @_;
#
#   my $padding = 0;  # 1=VCF, 0=RCF or Kaviar
#
#   # Initialize Kaviar if needed
#   if (ref($kav) ne 'Kaviar') {
#     $kav = new Kaviar($frz, $basedir);
#   }
#   exit if (ref($kav) ne 'Kaviar');
#
#   my $na = scalar @alleles;
#
#   # Require: no empty alleles
#   for (my $i = 0; $i < $na; $i++) {
#     return unless length($alleles[$i]);
#   }
#
#   my @split_alleles;
#   for (my $i = 0; $i < $na; $i++) {
#     my @chars = split //, $alleles[$i];
#     $split_alleles[$i] = \@chars;
#   }
#   my $trimmed = 0;
#
#   # while change in alleles do
#   my $change;
#   do {
#     $change = 0;
#     # if alleles end with same nucleotide then
#     #   truncate rightmost nucleotide of each allele
#     my $right_same = 1;
#     my $j = scalar @{$split_alleles[0]} - 1;
#     my $rightchar = $split_alleles[0]->[$j];
#     for (my $i = 1; $i < $na; $i++) {
#       my $j = scalar @{$split_alleles[$i]} - 1;
#       $right_same = 0 unless $split_alleles[$i]->[$j] eq $rightchar;
#     }
#     if ($right_same) {
#       for (my $i = 0; $i < $na; $i++) {
#         pop @{$split_alleles[$i]};
#       }
#       $change = 1;
#       $trimmed =  1;
#
#       # if there exists an empty allele then
#       #   extend alleles 1 nucleotide to left
#       my $any_empty = 0;
#       for (my $i = 0; $i < $na; $i++) {
#         $any_empty = 1 unless scalar @{$split_alleles[$i]};
#       }
#       if ($any_empty) {
#         $start--;
#         my $nt = $kav->refbase($chrom, $start);
#         for (my $i = 0; $i < $na; $i++) {
#           unshift @{$split_alleles[$i]}, $nt;
#         }
#       }
#     }
#   } until (!$change);
#
#   # while leftmost nucleotide of each allele are the same, cleave
#   my $left_same = 1;
#   while (1) {
#     # check whether we should cleave
#     last unless scalar @{$split_alleles[0]} > $padding;
#     my $leftchar = $split_alleles[0]->[0];
#     my $length_ok=1;
#     for (my $i = 1; $i < $na; $i++) {
#       $left_same = 0 unless (defined $split_alleles[$i]->[0]) &&
#           ($split_alleles[$i]->[0] eq $leftchar);
#           #$length_ok = 0 unless scalar((@{$split_alleles[$i]}));
#           $length_ok = 0 unless scalar((@{$split_alleles[$i]}) > $padding);
#     }
#     last unless $left_same && $length_ok;
#
#     # cleave
#     for (my $i = 0; $i < $na; $i++) {
#       shift @{$split_alleles[$i]};
#     }
#     $start++;
#     $trimmed = 1;
#   }
#
#   my @out_alleles;
#   for (my $i = 0; $i < $na; $i++) {
#     push @out_alleles, join('', @{$split_alleles[$i]});
#   }
#
#   return $start, @out_alleles if $trimmed;
#   return;
# }