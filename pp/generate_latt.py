# ==================
# GENERATE LATT FILE
# ==================

# python3 generate_latt.py BiVO4/BiVO4.xyz BiVO4/BiVO4-rates.dat 4.0
# python3 generate_latt.py Ga2O3/Ga2O3.xyz Ga2O3/Ga2O3-rates.dat 3.4

#Libraries
import numpy as np
import sys
    
# convert fractional to cartesian coordinates    
def frac2cart(v, p):
    a = p[0]; b = p[1]; c = p[2]; alpha = p[3]*np.pi/180;  beta = p[4]*np.pi/180; gamma = p[5]*np.pi/180;
    omega = a*b*c*np.sqrt( 1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    Mfc=[[a,    b*np.cos(gamma),  c*np.cos(beta)                                                ],
         [0.0,  b*np.sin(gamma),  c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/(np.sin(gamma))  ], 
         [0.0,  0.0,              omega/(a*b*np.sin(gamma))                                     ]]
    return np.dot(Mfc, v)

class Atom:
    def __init__(self, kind, r):
        """
        kind:  name atom
        coord: [x, y, z]
        """
        self.kind  = kind
        self.r = r
        self.NN_kind = []
        self.NN_cell = []
        self.NN_dist = []

class Rate:
    def __init__(self,trans,d,Ea,prefac):
        """
        trans: label transition involved (removing numbers)
        d:     distance between sites
        Ea:    activation energy
        """
        self.trans = ''.join([k for k in trans if not k.isdigit()])
        self.d = d
        self.Ea = Ea
        self.prefac = prefac

class Structure:
    """
    atoms: list of Atoms
    """
    # Initialization
    def __init__ (self, infile, infile_rates):

        # Parse structure
        with open(infile) as f:
            data = f.readlines()
        self.p = [float(x) for x in data[1].split()]
        self.atoms = []
        for line in data[2:]:
            kind = line.split()[0]
            r = [float(x) for x in line.split()[1:]]
            self.atoms.append(Atom(kind, r))
        self.N = len(self.atoms)
        self.a_vec = frac2cart([1,0,0],self.p)
        self.b_vec = frac2cart([0,1,0],self.p)
        self.c_vec = frac2cart([0,0,1],self.p)

        # mapping
        self.idx_atoms = {}
        for i in range(self.N):
            self.idx_atoms[self.atoms[i].kind] = i

        # Parse rates
        self.rates = []
        with open(infile_rates) as f:
            data = f.readlines()
        for line in data[1:]:
            trans = line.split()[0]
            [d,Ea,prefac] = [float(x) for x in line.split()[1:]]
            self.rates.append(Rate(trans,d,Ea,prefac))

    def dist(self,i,j,Nr=1):
        """
        calculate shortest distance between atom i and atom j
        accounting for Nr replicas on all directions
        """
        r1 = self.atoms[i].r
        r2 = self.atoms[j].r
        vals = [int(x) for x in np.linspace(-Nr,Nr,2*Nr+1)]
        dist = {} 
        for l in vals:
            for m in vals:
                for n in vals:
                    r2r = list(np.array(r2)+ l*np.array(self.a_vec) + m*np.array(self.b_vec) + n*np.array(self.c_vec))
                    key = "{:2d} {:2d} {:2d}".format(l,m,n)
                    dist[key] = np.sqrt(np.sum([np.abs(r1[k]-r2r[k])**2 for k in range(3)]))
        
        idx_min_d = np.argmin([dist[k] for k in dist.keys()])
        idx_cell = list(dist.keys())[idx_min_d]
        dmin = dist[idx_cell]
        return [dmin, idx_cell]

    def nearest_neighbour(self, dcut = float(sys.argv[3])):
        """
        return list of nearest neighbour within a cutoff distance dcut
        and return a list of indexes specifying the supercell repetition
        """
        for i in range(self.N):
            NN_kind = []; NN_cell = []; NN_dist = []
            for j in range(self.N):
                if j != i:
                    [d,idx_cell] = self.dist(i,j)
                    if d <= dcut:
                        NN_kind.append(self.atoms[j].kind)
                        NN_cell.append(idx_cell)
                        NN_dist.append(d)
            self.atoms[i].NN_kind = NN_kind
            self.atoms[i].NN_cell = NN_cell
            self.atoms[i].NN_dist = NN_dist

    def print_latt(self):
        """
        print input file for the random walk simulation
        """
        with open(sys.argv[1][:-4]+".latt", 'w') as f:
            f.write("3D\n")
            f.write("{:12.8f}  {:12.8f}  {:12.8f} \n".format(self.a_vec[0],self.a_vec[1],self.a_vec[2]))
            f.write("{:12.8f}  {:12.8f}  {:12.8f} \n".format(self.b_vec[0],self.b_vec[1],self.b_vec[2]))
            f.write("{:12.8f}  {:12.8f}  {:12.8f} \n".format(self.c_vec[0],self.c_vec[1],self.c_vec[2]))
            f.write("{:2d} {:2d} \n".format(self.N,1))
            for i in range(self.N):
                A = self.atoms[i]
                f.write("{:7s} {:4s} {:12.8f} {:12.8f} {:12.8f} {:4d} {:3s}".format(
                          A.kind, "N"+str(i+1), A.r[0],A.r[1],A.r[2],len(A.NN_kind), " "))
                for j in range(len(A.NN_kind)):
                     f.write("{:7s} {:4s} {:10s} {:2s}".format(A.NN_kind[j], "P1", A.NN_cell[j], " "))
                f.write("\n")
            f.write("864 \n")
            N_found = 0; N_not_found = 0; missing_hoppings = [] 
            for i in range(self.N):
                A = self.atoms[i]
                for j in range(len(A.NN_kind)):

                    # find the right transition energies for (i,inn) by selecting the right type and distance
                    type_i = ''.join([k for k in A.kind       if not k.isdigit()]) 
                    type_j = ''.join([k for k in A.NN_kind[j] if not k.isdigit()]) 
                    trans = type_i+"-"+type_j
                    Ea = 0; prefac = 0
                    
                    for R in self.rates:
                        if R.trans == trans and abs(R.d - A.NN_dist[j]) < 0.05:
                            Ea = R.Ea
                            prefac = R.prefac
                            break
                    if (Ea == 0 and prefac == 0):
                        N_not_found += 1
                        missing_hoppings.append("{:10s} d = {:4.2f}".format(type_i+"-"+type_j,round(A.NN_dist[j],4)))
                    else:
                        N_found += 1

                    # symmetrize to find the opposite transition
                    Ea2 = 0
                    trans2 = type_j+"-"+type_i
                    for R in self.rates:
                        if R.trans == trans2 and abs(R.d - A.NN_dist[j]) < 0.01:
                            Ea2 = R.Ea
                            break
                    
                    f.write("{:4s} {:4s} {:4s} {:6.4f} {:6.4f} {:6.4f} \n".format(
                            "N"+str(i+1), "N"+str(self.idx_atoms[A.NN_kind[j]]+1), "P1", Ea, Ea2, prefac))
            print("{:d} transitions: {:d} found, {:d} not found".format(N_found+N_not_found,N_found, N_not_found))
            
            # Print missing hoppings
            missing_hoppings = list(dict.fromkeys(missing_hoppings))
            if len(missing_hoppings) > 0: 
                print("Missing hoppings:")
                for s in missing_hoppings:
                    print(s)
    
# Inputs
S = Structure(sys.argv[1],sys.argv[2])
S.nearest_neighbour()
S.print_latt()