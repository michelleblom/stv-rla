import sys
import os

if __name__ == "__main__":

    with open(sys.argv[1], 'r') as listf:
        lines = [l.strip() for l in listf.readlines()]

        for f in lines:
            print(f)
            with open(f, 'r') as data:
                flines = data.readlines()


                q = flines[2].strip()
                v = flines[3].strip()

                o1 = flines[4].strip()
                o2 = flines[5].strip()

                path,name = os.path.split(f)
               
                newname = name
                if ".blt.csv" in name: 
                    newname = name.replace(".blt.csv", "")
                else:    
                    newname = name.replace(".csv", "")

                newname = newname.replace("summary_", "")

                otfile = path + "/" + newname + ".otc"
                print(otfile)
                with open(otfile, 'w') as ot:
                    print(o1, file=ot)
                    print(o2, file=ot) 

                qfile = path + "/" + newname + ".quota"
                print(qfile)
                with open(qfile, 'w') as qf:
                    print(q, file=qf)

               
                vfile = path + "/" + newname + ".voters"
                print(vfile)
                with open(vfile, 'w') as vf:
                    print(v, file=vf) 


