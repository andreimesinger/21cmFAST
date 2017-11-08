#!/usr/bin/perl


##### USAGE: extract_delTps.pl k_mode <delta_T file 1> <delta_T file 2> > outfile

my $pivot_k = shift(@ARGV);


sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

my $lpivot_k = log10($pivot_k);
my $ct=0;
my @lines;
foreach $file (@ARGV){
    $file =~ /\_nf(\d{1})\.(\d{6})\_/;
    $nf = "$1.$2";
    if ($nf < 0.01){
	next;
    }
    $file =~ /\_aveTb(\S+)\_(\S+)_(\S+)_(\S+)/;
    my $aveTb = $1;
    my $aveTbsq = $aveTb*$aveTb;
    $file =~ /\_z(\d{3})\.(\d{2})\_/;
    $z = "$1.$2";
    open IN, "<$file";
    while (<IN>){
	chomp;
	my ($k, $p, $perr) = split /\s/;
	my $lk = log10($k);
#	my $lp = log10($p/$aveTbsq); # dimentionless power
	my $lp = log10($p);
	if ($k > $pivot_k){
	    # get slope
	    my $dlk = $lk - $prev_lk;
	    my $dlp = $lp - $prev_lp;
	    my $slope = $dlp/$dlk;

	    # get interpolated power
	    my $lpivot_p = $prev_lp + $slope * ($lpivot_k - $prev_lk);
	    $lines[$ct++] = "$z\t$nf\t$aveTb\t$lpivot_k\t$lpivot_p\t$slope\n";
	    close IN;
	    last;
	}
	$prev_lk = $lk;
	$prev_lp = $lp;
 
   }
}

print sort(@lines);
