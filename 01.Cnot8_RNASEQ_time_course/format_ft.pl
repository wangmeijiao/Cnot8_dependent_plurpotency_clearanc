use strict;
use warnings;



while(<stdin>){
  chomp;
  next if($_=~/^#/) ;
  my @box = split/[\t ]+/;
  my @chrs = split/;/,$box[1];
  if(&check_same(\@chrs) == 0 ){ print STDERR "chr diffs\t$_\n";next  }
  my @starts=split/;/,$box[2];
  #if(&check_increase(\@starts) == 0 ){ print STDERR "start not increase\t$_\n";next  }
  my @ends = split/;/,$box[3];
  #if(&check_increase(\@ends) == 0 ){ print STDERR "end not increase\t$_\n";next  }
  #my ($flag, $sum_len) = &check_length(\@starts, \@ends, $box[5]);
  #if($flag == 0 ){ print STDERR "length diff($sum_len:$box[5])\t$_\n";next  }
  my @strands=split/;/,$box[4];
  if(&check_same(\@strands) == 0 ){ print STDERR "strand diffs\t$_\n";next  }
  my $data = join "\t", @box[5..$#box];
  print "$box[0]\t$chrs[0]\t$starts[0]\t$ends[-1]\t$data\n";

}

##sub

sub check_same(){
  my $idx = shift;
  my $chr;
  my $chr_flag = 1;
  foreach my $ctrl(@$idx){
     if(!defined $chr){ $chr = $ctrl}else{
       if($chr ne $ctrl ){$chr_flag=0;last}
     }
  }
  return $chr_flag;
}

sub check_increase(){
  my $idx = shift;
  my $chr;
  my $chr_flag = 1;
  foreach my $ctrl(@$idx){
     if(!defined $chr){ $chr = $ctrl}else{
       if($chr >= $ctrl ){$chr_flag=0;last}
     }
  }
  return $chr_flag;

}

sub check_length(){
   my ($starts,$ends,$length) = @_;
   #my @starts = split/,/,$starts;
   #my @ends = split/,/,$ends;
   die  if(scalar @$starts != scalar @$ends);
   my $sum;
   for (0..$#$starts){ $sum+= $ends->[$_] - $starts->[$_] }
   if($sum == $length){return (1,$sum )  }else{ return (0, $sum) }   
}



