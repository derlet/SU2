# This is a copy of the template file '/home/derlet/Software/ovito-pro-3.2.0-x86_64/share/ovito/scripts/modifiers/SU2.py'.
# Feel free to modify the code below as needed.

# This is a copy of the template file '/home/derlet/Software/ovito-pro-3.2.0-x86_64/share/ovito/scripts/modifiers/SU2.py'.
# Feel free to modify the code below as needed.

#
# SU(2) analysis Ovito modifyer, version 1.0, by Peter Derlet, see reference:
#
# P. M. Derlet, Correlated disorder in a well relaxed model binary glass 
# through a local SU(2) bonding topology, arXiv:2007.08878 (2020).
#
# The author acknowledges help given by Alexander Stukowski and Constanze Kalcher at www.ovito.org
#
from ovito.data import *
from ovito.modifiers import *
import numpy as np
#
# This helper function takes a two-dimensional array and computes a frequency
# histogram of the data rows using some NumPy magic.useRadii
# It returns two arrays (of equal length):
#    1. The list of unique data rows from the input array
#    2. The number of occurences of each unique row
# Adapted from OVITO coding example by Alexander Stukowski (See OVITO documentation on python coding)
#
def row_histogram(a):
    ca = np.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = np.unique(ca, return_index=True, return_inverse=True)
    counts = np.bincount(inverse)
    return (a[indices],counts)
#
# This user-defined modifier function gets automatically called by OVITO whenever the data pipeline is newly computed.
# It receives two arguments from the pipeline system:
# 
#    frame - The current animation frame number at which the pipeline is being evaluated.
#    data   - The DataCollection passed in from the pipeline system. 
#                The function may modify the data stored in this DataCollection as needed.
# 
def modify(frame, data,radicalVoronoi=True,modifiedVoronoi=True,maxCount=5):
    data.apply(VoronoiAnalysisModifier(use_radii = radicalVoronoi, generate_bonds = True))
    #
    # Obtain positions and bonding data
    #
    positions = data.particles['Position']
    bond_topology = data.particles.bonds.topology  
    bonds_enum = BondsEnumerator(data.particles.bonds)
    #
    # Construct list of neighbours according Voronoi generated bonds
    #
    lnl=[]
    for particle_index in range(data.particles.count):
        yield(particle_index/data.particles.count)
        nl=[]
        for bond_index in bonds_enum.bonds_of_particle(particle_index):
            a = bond_topology[bond_index, 0]
            b = bond_topology[bond_index, 1]
            if a==particle_index:
                nl.append(b)
            else:
                nl.append(a)
        lnl.append(nl)
    #
    # Statistics of bond order for standard Voronoi tessellation.
    # Calculates the global bond order and average coordination, and constructs a historgram of bond-orders.
    #
    hist=[]
    aq=0.0
    az=0.0
    bondNumber=0
    for i_particle in range(data.particles.count):
        yield(i_particle/data.particles.count)
        az+=float(len(lnl[i_particle]))
        for j_particle in lnl[i_particle]:
            if j_particle<i_particle: continue
            bondOrder=len(list(set(lnl[i_particle]).intersection(set(lnl[j_particle]))))
            hist.append(np.array([bondOrder]))
            aq+=float(bondOrder)
            bondNumber+=1
    aq/=float(bondNumber)
    az/=float(data.particles.count)
    eq=12.0/(6.0-aq) #Prediction of Euler's theorem
    print('Common neighbour average bond order     %f' % aq)
    print('Voronoi determined average coordination %f' % az)
    print('Eueler theorem coordination             %f' % eq)
    #
    # Generate historgram of bond orders
    #
    unique_indices, counts = row_histogram(np.array(hist))  
    print ('Bond order histogram:')
    for i in range(len(unique_indices)):
        if (counts[i]>0): print("%i\t%i\t%.4f %%" % (unique_indices[i][0],counts[i],100.0*float(counts[i])/len(hist))) 
    #
    # Modified Voronoi tessellation: exlude bonds which pass through a triplet of neighbouring atoms
    #
    if modifiedVoronoi:
        excluded_bond_number=0
        count=0
        while (count==0 or (count>0 and excluded_bond_number>0 and count<maxCount)):
            excluded_bond_number=0   
            for i_particle in range(data.particles.count):
                yield(i_particle/data.particles.count)
                for j_particle in lnl[i_particle]:
                    if j_particle<i_particle: continue
                    #
                    # Construct list of common neighbours between particles 'i_particle' and 'j_particle'
                    #
                    cnl=list(set(lnl[i_particle]).intersection(set(lnl[j_particle])))
                    #    
                    # Cycle through all triplets of these common neighbours and remove the bond if it intersects 
                    # the plane of the triangle defined by the current triplet. 
                    # 
                    break_triplet_loop=False    
                    for i_triplet in cnl:
                        for j_triplet in cnl:
                            if j_triplet<=i_triplet: continue
                            for k_triplet in cnl:
                                if k_triplet<=j_triplet: continue 
                                ij_member=j_triplet in lnl[i_triplet]
                                ik_member=k_triplet in lnl[i_triplet]
                                jk_member=k_triplet in lnl[j_triplet]
                                if ij_member and ik_member and jk_member:
                                                                                                                
                                    vecb=data.cell.delta_vector(positions[j_triplet],positions[i_triplet])
                                    vecc=data.cell.delta_vector(positions[i_particle],positions[i_triplet])
                                    vecd=data.cell.delta_vector(positions[j_particle],positions[i_triplet])
                                    tv_abde=np.dot(vecb,np.cross(vecc,vecd))/6.0
    
                                    vecb=data.cell.delta_vector(positions[k_triplet],positions[j_triplet])
                                    vecc=data.cell.delta_vector(positions[i_particle],positions[j_triplet])
                                    vecd=data.cell.delta_vector(positions[j_particle],positions[j_triplet])
                                    tv_bcde=np.dot(vecb,np.cross(vecc,vecd))/6.0

                                    vecb=data.cell.delta_vector(positions[i_triplet],positions[k_triplet])
                                    vecc=data.cell.delta_vector(positions[i_particle],positions[k_triplet])
                                    vecd=data.cell.delta_vector(positions[j_particle],positions[k_triplet])
                                    tv_cade=np.dot(vecb,np.cross(vecc,vecd))/6.0

                                    if (tv_abde<0.0 and tv_bcde<0.0 and tv_cade<0.0) or (tv_abde>0.0 and tv_bcde>0.0 and tv_cade>0.0):
                                        excluded_bond_number+=1
                                        lnl[i_particle].remove(j_particle)
                                        lnl[j_particle].remove(i_particle)
                                        break_triplet_loop=True
                                if break_triplet_loop: break    
                            if break_triplet_loop: break
                        if break_triplet_loop: break
            print('Iteration %i, number of bonds removed: %i ' % (count+1,excluded_bond_number))
            count+=1
        #
        # Statistics of bond order for standard Voronoi tessellation.
        # Calculates the global bond order and average coordination, and constructs a historgram of bond-orders.
        #
        hist=[]
        aq=0.0
        az=0.0
        bondNumber=0
        for i_particle in range(data.particles.count):
            yield(i_particle/data.particles.count)
            az+=float(len(lnl[i_particle]))
            for j_particle in lnl[i_particle]:
                if j_particle<i_particle: continue
                bondOrder=len(list(set(lnl[i_particle]).intersection(set(lnl[j_particle]))))
                hist.append(np.array([bondOrder]))
                aq+=float(bondOrder)
                bondNumber+=1
        aq/=float(bondNumber)
        az/=float(data.particles.count)
        eq=12.0/(6.0-aq) #Prediction of Euler's theorem
        print('Common neighbour average bond order     %f' % aq)
        print('Voronoi determined average coordination %f' % az)
        print('Eueler theorem coordination             %f' % eq)
        #
        # Generate historgram of bond orders
        # 
        print ('Bond order histogram:')
        unique_indices, counts = row_histogram(np.array(hist))   
        for i in range(len(unique_indices)):
            if (counts[i]>0): print("%i\t%i\t%.4f %%" % (unique_indices[i][0],counts[i],100.0*float(counts[i])/len(hist)))
    #
    # Create bond order property of bonds (removed bonds are assigned a bond order of zero)
    #                  
    bo = data.bonds_.create_property('Bond Order', dtype=int, components=1)
    maxBondOrder=0
    with bo:
        for bond_index in range(data.particles.bonds.count):
            yield(bond_index/data.particles.bonds.count)
            a = bond_topology[bond_index, 0]
            b = bond_topology[bond_index, 1]
            if a in lnl[b] and b in lnl[a]:
                bondOrder=len(list(set(lnl[a]).intersection(set(lnl[b]))))
            else:
                bondOrder=0
            bo[bond_index]=bondOrder
            if bondOrder>maxBondOrder: maxBondOrder=bondOrder
    #
    # Determine local SU(2) topology, identifying and labelling those used in arXiv:2007.08878
    #
    label=[]
    label.append([0,12,0])
    label.append([0,12,2])
    label.append([0,12,3])
    label.append([0,12,4])
    label.append([0,12,5])
    label.append([0,12,6])
    label.append([1,10,2])
    label.append([1,10,3])
    label.append([1,10,4])
    label.append([1,10,5])
    label.append([1,10,6])
    label.append([1,10,7])
    label.append([2,8,0])
    label.append([2,8,1])
    label.append([2,8,2])
    label.append([2,8,3])
    label.append([2,8,4])
    label.append([2,8,5])
    label.append([2,8,6])
    label.append([2,8,7])
    label.append([2,8,8])
    label.append([3,6,0])
    label.append([3,6,1])
    label.append([3,6,2])
    label.append([3,6,3])
    label.append([3,6,4])
    label.append([3,6,5])
    label.append([3,6,6])
    label.append([3,6,7])
    label.append([3,6,8])
    label.append([3,6,9])
    #
    # The SU2 poperty of an atom ranges between 0 and 31. Here 0 to 30 corresponds to the above labelling scheme 
    # and 31 represents an unlabelled environment referred to as "other"
    #
    ll = data.particles_.create_property('SU2 Label', dtype=int, components=1)
    with ll:    
        for i_particle in range(data.particles.count):
            hist=[0]*(maxBondOrder+1)
            for j_particle in lnl[i_particle]:
                hist[len(list(set(lnl[i_particle]).intersection(set(lnl[j_particle]))))]+=1
            if (hist[4]+hist[5]+hist[6]==len(lnl[i_particle])) and (12-hist[4]+hist[6]==len(lnl[i_particle])):
                found_label=False
                for i_label in range(len(label)):
                    if (hist[4]==label[i_label][0]) and (hist[5]==label[i_label][1]) and (hist[6]==label[i_label][2]):
                        ll[i_particle]=i_label
                        found_label=True
                        break
                if not(found_label): ll[i_particle]=31
            else :
                ll[i_particle]=31
    #
    # Determine mean bond order for each atom, defining the corresponding atom property
    #
    q = data.particles_.create_property('mean bond order', dtype=float, components=1)
    with ll:    
        for i_particle in range(data.particles.count):
            hist=[0]*(maxBondOrder+1)
            for j_particle in lnl[i_particle]:
                hist[len(list(set(lnl[i_particle]).intersection(set(lnl[j_particle]))))]+=1
            num=0.0
            den=0.0 
            for i in range(len(hist)):
                num+=hist[i]*i
                den+=hist[i]
            q[i_particle]=float(num)/float(den)
    #
    # Colour bonds according to bond order, where:
    # white  - removed bond
    # purple - 1,2,3 fold bond
    # green  - 4 fold bond
    # blue   - 5 fold bond
    # red    - 6 fold bond
    # grey   - 7+ fold bond 
    #
    data.bonds_.create_property('Color')
    color_scheme = ((1.0,1.0,1.0),(0.8,0.6,1.0),(0.8,0.6,1.0),(0.5,0,0.5),(0,1.0,0),(0,0,1.0),(1.0,0,0),(0.69,0.77,0.89),(0.69,0.77,0.89))
    for value in np.unique(bo):
        if value < len(color_scheme):
            data.bonds_['Color'][bo == value] = [color_scheme[value]]
        else:
            data.bonds_['Color'][bo == len(color_scheme)-1] = [color_scheme[len(color_scheme)-1]]


      
          