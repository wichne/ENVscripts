#!/usr/bin/perl

use lib $ENV{SCRIPTS};
use ENV;
use strict;
use Getopt::Std;
use Spreadsheet::ParseXLSX;
use Spreadsheet::ParseExcel;
use Excel::Writer::XLSX;

my %arg;
&getopts('f:D:p:u:o:',\%arg);
my $user = $arg{'u'} ? $arg{'u'} : $ENV{USER};
my $outfile = $arg{'o'} ? $arg{'o'} : $arg{'f'} . ".translated.xlsx";
my $xlsxout = Excel::Writer::XLSX->new($outfile);

my $dbh = connect(\%arg);

$| = 1;
my $parser = Spreadsheet::ParseXLSX->new;
my $wkbk = $parser->parse($arg{f});

# go through each worksheet in the workbook
for my $wksht ($wkbk->worksheets() ) {
    print STDERR "Found worksheet ", $wksht->get_name, "\n";

    # make a worksheet in the output workbook
    my $outSheet = $xlsxout->add_worksheet($wksht->get_name);

    # Now step through the rows
    my ( $row_min, $row_max ) = $wksht->row_range();
    my ( $col_min, $col_max ) = $wksht->col_range();

    # the first row is headers, or should be.
    # Look for text suggesting the column holds accessions (rather than annotation)
    my @id_cols;
    for my $col ($col_min .. $col_max) { 
	my $cell = $wksht->get_cell($row_min, $col);
	if (! defined $cell) { print STDERR "No cell: $row_min, $col\n"; next; }
	if ($cell->value =~ /[Bb]in|HL/) {
	    push @id_cols, $col;
	}
    }
    for my $row ( $row_min+1 .. $row_max ) {
	my $fig_prefix;
	for my $col ( @id_cols) {
	    my $cell = $wksht->get_cell($row, $col);
	    if (! defined $cell) {
#		print STDERR "No cell: $row_min, $col\n";
		next;
	    }
	    my $value = $cell->value;
	    # my $format = translate_format($cell->get_format, $xlsxout);
	    my $format;

	    $value =~ s/^\s*//;
	    $value =~ s/\s*$//;
	    my @ids = split/[\s\,]+/, $value; 
	    my @locus;
	    my $format = $xlsxout->add_format;
	    foreach my $i (@ids) {
		$i =~ s/fig\|//;
		my $fid = get_feature_id_by_accession($dbh, $i);
		my $locusr = get_locus_tag_by_feature_id($dbh, $fid);
		if (!@$locusr) {
		    print STDERR "!!!!PROBLEM: Couldn't find locus for $fid $i.\n";
		    $format->set_bg_color('red');
		    push @locus, $i;
		} else {
		    push @locus, @$locusr;
		}
	    }
	    my $new_value = join(",", @locus);
	    $outSheet->write($row, $col, $new_value, $format);
	}
    }
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
