# VMD/Tcl: extract molecules within 6 Å and 20 Å of each albumin.
# Requires an albumin segment list, one segment per line, e.g. inputs/albumin_segs.txt.

set psf_file "/path/to/dry.psf"
set dcd_template "/path/to/dry.%d.dcd"
set dcd_start 49
set dcd_end 63
set dcd_stride 25
set output_dir "."
set albumin_seg_file "inputs/albumin_segs.txt"
set cutoffs {6 20}

mol new $psf_file type psf waitfor all
for {set i $dcd_start} {$i <= $dcd_end} {incr i} {
    mol addfile [format $dcd_template $i] type dcd step $dcd_stride waitfor all
}
set num_frames [molinfo top get numframes]

set albumin_segs {}
set fh [open $albumin_seg_file r]
while {[gets $fh line] >= 0} {
    set line [string trim $line]
    if {$line ne ""} { lappend albumin_segs $line }
}
close $fh

foreach albumin_seg $albumin_segs {
    set albumin_dir [file join $output_dir $albumin_seg]
    if {![file exists $albumin_dir]} { file mkdir $albumin_dir }
    set target_selection "segid $albumin_seg"

    for {set frame 0} {$frame < $num_frames} {incr frame} {
        animate goto $frame
        foreach cutoff $cutoffs {
            set sel [atomselect top "same residue as within $cutoff of ($target_selection)"]
            if {[llength [$sel get index]] > 0} {
                if {$cutoff == 6} {
                    set filename [format "Albumin_%s_frame_%04d.pdb" $albumin_seg $frame]
                } else {
                    set filename [format "Albumin_%s_%dA_frame_%04d.pdb" $albumin_seg $cutoff $frame]
                }
                set filepath [file join $albumin_dir $filename]
                $sel writepdb $filepath
                puts "Wrote $filepath"
            }
            $sel delete
        }
    }
    puts "Completed albumin $albumin_seg"
}

# quit
