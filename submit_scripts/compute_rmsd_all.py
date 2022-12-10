import sys
sys.path.insert(0,'/home/users/jwang003/docking')
import compute_rmsds_all

compute_rmsds_all.struct_align_all(sys.argv[1],sys.argv[2:])
print('aligned')
compute_rmsds_all.struct_sort(sys.argv[1:])
print('split')
compute_rmsds_all.rmsd_all()
