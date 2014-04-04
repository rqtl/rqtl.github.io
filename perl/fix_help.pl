#!/usr/bin/perl
######################################################################
# fix the html files for the R/qtl help on its web site
#   -change the name of est.map.html to est_map.html
#   -fix links in all files
######################################################################

$dir = "manual/html";
@oldfile = ("est.map.html", "plot.map.html", "pull.map.html", 
	    "replace.map.html", "sim.map.html", "summary.map.html");
@newfile = ("est_map.html", "plot_map.html", "pull_map.html", 
	    "replace_map.html", "sim_map.html", "summary_map.html");

$n = @oldfile - 1;
foreach $i (0..$n) {
    if(-e "$dir/$oldfile[$i]") {
	system("\\mv $dir/$oldfile[$i] $dir/$newfile[$i]");
	print("Moving $oldfile[$i] to $newfile[$i]\n");
    }
}
 
opendir(DIR, $dir);
while(defined($file = readdir(DIR))) {
    unless($file =~ /.html$/) {
	next; 
    }
    $ifile = "$dir/$file";
    open(IN, $ifile) or next;
    $ofile = "$dir/junk";
    open(OUT, ">$ofile") or die("Cannot write to junk");
    $flag = 0;
    while($line = <IN>) {
	foreach $i (0..$n) {
	    if($line =~ /$oldfile[$i]/) {
		$line =~ s/$oldfile[$i]/$newfile[$i]/g;
		$flag = 1;
	    }
	    if($line =~ /Sen, \'S\./) {
		$line =~ s/Sen, \'S./Sen, &\#346;./;
		$flag = 1;
	    }
	}
	print OUT $line;
    }
    close(OUT);
    close(IN);
    if($flag) {
	system("\\cp $ofile $ifile");
	print("Copying $ofile to $ifile\n");
    }
}
if(-e $ofile) {
    system("\\rm $ofile");
}

    
