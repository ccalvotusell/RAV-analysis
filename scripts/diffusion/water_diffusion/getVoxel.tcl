animate delete beg 0 end 0 waitfor all
animate goto end
sleep 1
set output "/.../Analysis/Mass_voxel/pdbs"
for {set x 200} {$x <= 2800} {incr x 200} {
	for {set y 200} {$y <= 2800} {incr y 200} {
		for {set z 200} {$z <= 2800} {incr z 200} {
			set sel [atomselect top "x < $x + 100 and x > $x - 100 and y < $y + 100 and y > $y - 100 and z < $z + 100 and z > $z - 100"]
			set num [$sel num]
			puts "$x $y $z"
			puts "$num"
			if {$num > 10000} {
				$sel writepdb ${output}/x${x}_y${y}_z${z}.pdb
				$sel writepsf ${output}/x${x}_y${y}_z${z}.psf
				animate write dcd ${output}/x${x}_y${y}_z${z}.dcd beg 0 end -1 waitfor all sel $sel 
			}
			$sel delete
		}
	}
}
			
