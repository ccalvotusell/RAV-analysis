# VMD/Tcl: extract molecules within 6 Å and 20 Å of each SARS-CoV-2 spike.
# Run from the desired output directory, e.g.:
#   vmd -dispdev text -e scripts/01_extract_spike_contact_shells.tcl

# -----------------------------
# User settings
# -----------------------------
set psf_file "/path/to/dry.psf"
set dcd_template "/path/to/dry.%d.dcd"
set dcd_start 49
set dcd_end 63
set dcd_stride 25
set output_dir "."
set spike_ids {}

# Default: spike00-spike29, excluding spike06.
for {set n 0} {$n <= 29} {incr n} {
    if {$n != 6} { lappend spike_ids [format "%02d" $n] }
}

set cutoffs {6 20}
set output_prefix "RAV_v27"

# -----------------------------
# Load trajectory
# -----------------------------
mol new $psf_file type psf waitfor all
for {set i $dcd_start} {$i <= $dcd_end} {incr i} {
    set dcd_file [format $dcd_template $i]
    mol addfile $dcd_file type dcd step $dcd_stride waitfor all
}
set num_frames [molinfo top get numframes]
puts "Loaded $num_frames frames."

# -----------------------------
# Extract contact shells
# -----------------------------
foreach spike_num $spike_ids {
    set target_selection "segid A1${spike_num} B1${spike_num} C1${spike_num} A2${spike_num} B2${spike_num} C2${spike_num} ${spike_num}*"
    set spike_dir [file join $output_dir [format "spike%s" $spike_num]]
    if {![file exists $spike_dir]} { file mkdir $spike_dir }

    for {set frame 0} {$frame < $num_frames} {incr frame} {
        animate goto $frame
        foreach cutoff $cutoffs {
            set sel [atomselect top "same residue as within $cutoff of ($target_selection)"]
            if {[llength [$sel get index]] > 0} {
                if {$cutoff == 6} {
                    set filename [format "%s_s%s_frame_%04d.pdb" $output_prefix $spike_num $frame]
                } else {
                    set filename [format "%s_s%s_%dA_frame_%04d.pdb" $output_prefix $spike_num $cutoff $frame]
                }
                set filepath [file join $spike_dir $filename]
                $sel writepdb $filepath
                puts "Wrote $filepath"
            } else {
                puts "No atoms within $cutoff Å for spike $spike_num at frame $frame"
            }
            $sel delete
        }
    }
    puts "Completed spike $spike_num"
}

# quit
