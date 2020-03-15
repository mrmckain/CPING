#!/usr/env perl


my %gff;
open my $gff, "<", $ARGV[0];
while(<$gff>){
		chomp;
		if(/\#\#/){
			next;
		}
		my @tarray = split /\s+/;

		if($tarray[2] eq "gene"){
			$gff{$tarray[0]}{$tarray[3]+$ARGV[2]}={$tarray[4]+$ARGV[2]};
		}
}

open my $seq, "<", $ARGV[1];

my %genome;
my $sid;
while(<$seq>){
		chomp;
		if(/>/){
			$_ =~ />(.*?)\s+.+$/;
			$sid=$1;
		}
		else{
			$genome{$sid}.=$_;
		}
}


my %alt_nuc;
push(@{$alt_nuc{K}}, "G");
push(@{$alt_nuc{K}}, "T");
push(@{$alt_nuc{M}}, "A");
push(@{$alt_nuc{M}}, "C");
push(@{$alt_nuc{R}}, "A");
push(@{$alt_nuc{R}}, "G");
push(@{$alt_nuc{Y}}, "C");
push(@{$alt_nuc{Y}}, "T");
push(@{$alt_nuc{S}}, "C");
push(@{$alt_nuc{S}}, "G");
push(@{$alt_nuc{W}}, "A");
push(@{$alt_nuc{W}}, "T");
push(@{$alt_nuc{B}}, "C");
push(@{$alt_nuc{B}}, "G");
push(@{$alt_nuc{B}}, "T");
push(@{$alt_nuc{V}}, "C");
push(@{$alt_nuc{V}}, "G");
push(@{$alt_nuc{V}}, "A");
push(@{$alt_nuc{H}}, "C");
push(@{$alt_nuc{H}}, "A");
push(@{$alt_nuc{H}}, "T");
push(@{$alt_nuc{D}}, "A");
push(@{$alt_nuc{D}}, "G");
push(@{$alt_nuc{D}}, "T");
push(@{$alt_nuc{D}}, "A");
push(@{$alt_nuc{D}}, "T");
push(@{$alt_nuc{D}}, "G");
push(@{$alt_nuc{D}}, "C");

#my @enz = split(/,/, $resenz);
open my $enz_files, "<", $ARGV[3];
my %enzymes;
my %used;
while(<$enz_files>){
	chomp;
	my @tarray = split/\s+/;
	if(exists $used{$tarray[2]}){
		next;
	}
	$enzymes{$tarray[1]}=1;
	$used{$tarray[2]}=1;

}

#my %enzymes = map {$_ => 1} @enz;
my %sites;

open my $enz_files, "<", $ARGV[3];
while(<$enz_files>){
		chomp;
		my @tarray = split /\s+/;
		if(exists $enzymes{$tarray[1]}){
			my @possible_sites = site_generator($tarray[0]);
			for my $psite (@possible_sites){
				push(@{$sites{$tarray[1]}}, $psite);
				my $tempsite = $psite;
				$tempsite =~ tr/ATCGatcg/TAGCtagc/;
				$tempsite = reverse($tempsite);
				push(@{$sites{$tarray[1]}}, $tempsite);
			}
		}
}

my %coding;
my %noncoding;
my $start_0=0;
for my $chr (keys %gff){
		for my $start (sort {$a<=>$b} keys %{$gff{$chr}}){
			$noncoding{$chr}{$start_0}{$start-2}=substr($genome{$chr}, $start_0, ($start-$start_0));
			$coding{$chr}{$start-1}{$gff{$chr}{$start}-1}=substr($genome{$chr}, $start-1, ($gff{$chr}{$start}-$start-1));
		}
}
			

my %positions;
for my $ids (sort keys %coding){
	for my $start (keys %{$coding{$ids}}){
		for my $end (keys %{$coding{$ids}{$start}}){
	
	my $offset = 0;
	for my $enzyme (sort keys %sites){
		for my $cut (@{$sites{$enzyme}}){
			$offset = 0;
			#print "cut is $cut\n";
			while($offset<length($coding{$ids}{$start}{$end})){
				my $pos = index($coding{$ids}{$start}{$end}, $cut, $offset);
			#	print "$pos\n";
				if($pos==-1){
					$offset= length($coding{$ids}{$start}{$end});
				}
				unless($pos==-1){
					$positions{$enzyme}{coding}++;
					$offset = $pos +1;
				}
			}
		}
	}
}}}

for my $ids (sort keys %noncoding){
	for my $start (keys %{$noncoding{$ids}}){
		for my $end (keys %{$noncoding{$ids}{$start}}){
	
	my $offset = 0;
	for my $enzyme (sort keys %sites){
		for my $cut (@{$sites{$enzyme}}){
			$offset = 0;
			#print "cut is $cut\n";
			while($offset<length($noncoding{$ids}{$start}{$end})){
				my $pos = index($noncoding{$ids}{$start}{$end}, $cut, $offset);
			#	print "$pos\n";
				if($pos==-1){
					$offset= length($noncoding{$ids}{$start}{$end});
				}
				unless($pos==-1){
					$positions{$enzyme}{noncoding}++;
					$offset = $pos +1;
				}
			}
		}
	}
}}}



my %enzyme_cuts;



my %redux_enz;

for my $enz1 (keys %positions){
	my $diff = $positions{$enz1}{coding}-$positions{$enz1}{noncoding};
	$redux_enz{$enz1}{$diff}=1;
}

my %top_three;
my $counter=0;
for my $enz1 (keys %redux_enz){
	for my $diff (sort {$b<=>$a} keys $redux_enz{$diff}){
		if($counter>=3){
			next;
		}
		$top_three{$enz1}=$diff;
		$counter++;
	}
}


open my $OUT, ">", "$ARGV[1]" . "_best_enzymes.txt";
print $OUT "Enzyme\tCoding_Cuts\tNoncoding_Cuts\n";
for my $enz1 (keys %top_three){

	print $OUT "$enz1\t$enzyme_cuts{$enz1}{coding}\t$enzyme_cuts{$enz1}{noncoding}\n";
}

###########################
sub site_generator {
	
	my $seq = $_[0];

	my %temp_seqs;

	my @split_seq = split(//, $seq);
	for my $base (@split_seq){
			if(exists $alt_nuc{$base}){
					if(%temp_seqs){
							for my $tseq (keys %temp_seqs){
								for my $nbase (@{$alt_nuc{$base}}){
									my $nseq .= $tseq . $nbase;
									$temp_seqs{$nseq}=1;
								}
								delete $temp_seqs{$tseq};
							}
					}
					else{
						for my $nbase (@{$alt_nuc{$base}}){
							$temp_seqs{$nbase}=1;
						}
					}
			}
			else{
				if(%temp_seqs){
						for my $tseq (keys %temp_seqs){
							my $nseq .= $tseq . $base;
							delete $temp_seqs{$tseq};
							$temp_seqs{$nseq}=1;
						}
				}
				else{
					$temp_seqs{$base}=1;
				}

			}
	}
	my @final_seqs;
	for my $seqs (keys %temp_seqs){
		push(@final_seqs, $seqs);
	}
	return @final_seqs;
}
