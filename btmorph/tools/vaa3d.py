"""
Script to convert Vaa3D SWC files to valid three-point soma files
"""

import sys,time

def convert_vaa3d(f_name,out_f):
    f = open(f_name,"r")
    orig_data = {}
    for line in f:
        if not line.startswith("#"):
            split = line.split()
            
            index = int(split[0].rstrip())
            swc_type = int(split[1].rstrip())
            x = float(split[2].rstrip())
            y = float(split[3].rstrip())
            z = float(split[4].rstrip())
            radius = float(split[5].rstrip())
            parent_index = int(split[6].rstrip())

            orig_data[index] = (swc_type,x,y,z,radius,parent_index)

    # Now write a new SWC description
    o = open(out_f,"w")
    t_w = ""
    # make the soma. vaa3d uses a fake one-point soma
    o_r = orig_data[1][4]
    t_w = t_w+"1 1 {0} {1} {2} {3} -1\n".format(orig_data[1][1],\
                                                        orig_data[1][2],\
                                                        orig_data[1][3],\
                                                        orig_data[1][4])
    t_w = t_w+"1 1 {0} {1} {2} {3} 1\n".format(orig_data[1][1],\
                                                        orig_data[1][2]-o_r,\
                                                        orig_data[1][3],\
                                                        orig_data[1][4])
    t_w = t_w+"1 1 {0} {1} {2} {3} 1\n".format(orig_data[1][1],\
                                                        orig_data[1][2]+o_r,\
                                                        orig_data[1][3],\
                                                        orig_data[1][4])
    o.write(t_w)

    # process all dendrites
    t_w=""
    mapping={}
    new_index=4
    for old_index in sorted(orig_data.keys())[1:]:
        #print "old_index: ", old_index
        mapping[old_index] = new_index
        t_w = t_w+"".format(new_index)

    # clean and leave
    o.flush()
    o.close()

if __name__=="__main__":
    f_name = sys.argv[1]
    out_f = f_name.split(".swc")[0]+"_clean.swc"
    print("Converting file {0} to {1}".format(f_name,out_f))
    convert_vaa3d(f_name,out_f)

