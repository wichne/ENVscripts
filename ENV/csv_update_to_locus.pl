#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;

my %arg;
&getopts('f:D:p:u:o:',\%arg);
my $user = $arg{'u'} ? $arg{'u'} : $ENV{USER};
my $outfile = $arg{'o'} ? $arg{'o'} : $arg{'f'} . ".translated";

my $dbh = connect(\%arg);

$| = 1;
my $parser;

open my $IN, $arg{'f'};
open(OUT, ">$arg{o}");

while (my $line = <$IN>) {
    chomp $line;
    my @acc = split/\t/, $line;    

    # the first row is headers, or should be.
    # Look for text suggesting the column holds accessions (rather than annotation)
    for (my $i=0; $i<@acc; $i++) {
	$acc[$i] =~ s/^\s*//;
	$acc[$i] =~ s/\s*$//;
	my @ids = split/[\s\,\|\;]+/, $acc[$i]; 
	my @locus;
	foreach my $i (@ids) {
	    if ($i =~ /fig/) { next }
	    $i =~ s/^gnl\_//;
	    my $fid = get_feature_id_by_accession($dbh, $i);
	    my $locusr = get_locus_tag_by_feature_id($dbh, $fid);
	    if (!@$locusr) {
		print STDERR "!!!!PROBLEM: Couldn't find locus for $fid $i.\n";
		push @locus, $i;
	    } else {
		push @locus, @$locusr;
	    }
	}
	my $new_value = join(",", @locus);
	print STDERR "$acc[$i] becomes $new_value\n";
	$acc[$i] = $new_value;
    }
    print OUT join("\t", @acc) . "\n";
}

sub translate_format {
    my $AlignHref = {0 => '',
		     1 => 'left',
		     2 => 'center',
		     3 => 'right',
		     4 => 'fill',
		     5 => 'justify',
		     6 => 'center_across',
		     7 => 'fill'};
    my $AlignVref = {0 => 'top',
		     1 => 'vcenter',
		     2 => 'bottom',
		     3 => 'vjustify',
		     4 => 'center'};
    my $old_format = shift;
    my $new_wkbk = shift;
    my $new_format = $new_wkbk->add_format;
    
    if (defined $old_format->{Font}) {
	my $fontobj = $old_format->{Font};
	if (defined $fontobj->{Name}) {
	    $new_format->set_font($fontobj->{Name});
	}
	if (defined $fontobj->{Bold}) {
	    $new_format->set_bold($fontobj->{Bold});
	}
	if (defined $fontobj->{Italic}) {
	    $new_format->set_italic($fontobj->{Italic});
	}
	if (defined $fontobj->{Height}) {
	    $new_format->set_size($fontobj->{Height});
	}
	if (defined $fontobj->{UnderlineStyle}) {
	    $new_format->set_underline($fontobj->{UnderlineStyle});
	} elsif (defined $fontobj->{Underline}) {
	    $new_format->set_underline($fontobj->{Underline});
	}
	if (defined $fontobj->{Color}) {
	    $new_format->set_color($fontobj->{Color});
	}
	if (defined $fontobj->{Strikeout}) {
	    $new_format->set_font_strikeout($fontobj->{Strikeout});
	}
	if (defined $fontobj->{Super}) {
	    $new_format->set_font_script($fontobj->{Super});
	}
    }
    if (defined $old_format->{AlignH}) {
	if ($old_format->{AlignH}) { $new_format->set_align($AlignHref->{$old_format->{AlignH}});}
    }
    if (defined $old_format->{AlignV}) {
	$new_format->set_align($AlignVref->{$old_format->{AlignV}});
    }
    if (defined $old_format->{Indent}) {
	$new_format->set_indent($old_format->{Indent});
    }
    if (defined $old_format->{Wrap}) {
	$new_format->set_text_wrap($old_format->{Wrap});
    }
    if (defined $old_format->{Shrink}) {
	$new_format->set_shrink($old_format->{Shrink});
    }
    if (defined $old_format->{Rotate}) {
	if ($old_format->{Rotate} == 1) {
	    $new_format->set_rotation(180);
	} elsif ($old_format->{Rotate} == 2) {
	    $new_format->set_rotation(-90);
	} elsif ($old_format->{Rotate} == 3) {
	    $new_format->set_rotation(90);
	} else { 
	    $new_format->set_rotation($old_format->{Rotate});
	}
    }
    if (defined $old_format->{JustLast}) {
	$new_format->set_text_justlast($old_format->{JustLast});
    }
#    if (defined $old_format->{ReadDir}) {
#    }
    if (defined $old_format->{BdrStyle}) {
	$new_format->set_left($old_format->{BdrStyle}->[0]);
	$new_format->set_right($old_format->{BdrStyle}->[1]);
	$new_format->set_top($old_format->{BdrStyle}->[2]);
	$new_format->set_bottom($old_format->{BdrStyle}->[3]);
    }
    if (defined $old_format->{BdrColor}) {
	$new_format->set_left_color($old_format->{BdrColor}->[0]);
	$new_format->set_right_color($old_format->{BdrColor}->[1]);
	$new_format->set_top_color($old_format->{BdrColor}->[2]);
	$new_format->set_bottom_color($old_format->{BdrColor}->[3]);
    }
    if (defined $old_format->{BdrDiag}) {
	$new_format->set_diag_type($old_format->{BdrDiag}->[0]);
	$new_format->set_diag_border($old_format->{BdrDiag}->[1]);
	$new_format->set_diag_color($old_format->{BdrDiag}->[2]);
    }
    if (defined $old_format->{Fill}) {
	$new_format->set_pattern($old_format->{Fill}->[0]);
	$new_format->set_fg_color($old_format->{Fill}->[1]);
	$new_format->set_bg_color($old_format->{Fill}->[2]);
    }
    if (defined $old_format->{Lock}) {
	$new_format->set_locked($old_format->{Lock});
    }
    if (defined $old_format->{Hidden}) {
	$new_format->set_hidden($old_format->{Hidden});
    }
    return $new_format;
}
