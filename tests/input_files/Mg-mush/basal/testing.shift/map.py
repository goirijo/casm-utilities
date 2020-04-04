import casmutils as cu
import glob

def main():
    files=glob.glob("shift*/cleave_0.00*/POSCAR")
    strucs=[cu.xtal._xtal.Structure.from_poscar(f) for f in files]
    s0=strucs[0]
    s1=strucs[0]

    print(s0)
    print(s1)

    cu.mapping._mapping.map_structure(s0,s1)

if __name__=="__main__":
    main()
