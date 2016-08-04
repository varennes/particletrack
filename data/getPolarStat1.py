import os
import numpy as np

# list all fort.1## files in current directory
onlyfiles = [f for f in os.listdir('.') if os.path.isfile(os.path.join('.', f))]
temp = []
for f in onlyfiles:
    if 'fort' in f:
        temp.append(f)
onlyfiles = temp
size1 = os.stat((onlyfiles[-1])).st_size
size2 = os.stat((onlyfiles[-2])).st_size
if size1 != size2:
    onlyfiles.pop()
print onlyfiles

pxM = []
pyM = []
pzM = []
pxV = []
pyV = []
pzV = []
pM  = []
pV  = []
e   = []
# loop through each fort.1## file
for f in onlyfiles:
    data = []
    # open and read fort.1## file
    with open( f, 'r') as fo:
        data = np.genfromtxt(fo);
        px = []
        py = []
        pz = []
        p  = []
        for x in data:
            px.append(x[0])
            py.append(x[1])
            pz.append(x[2])

            p.append( np.sqrt( x[0]**2 + x[1]**2 + x[2]**2) )

        pxM.append( np.mean(px))
        pyM.append( np.mean(py))
        pzM.append( np.mean(pz))
        pM.append( np.mean(p))

        pxV.append( np.std(px)**2)
        pyV.append( np.std(py)**2)
        pzV.append( np.std(pz)**2)
        pV.append( np.std(p)**2)

        e.append( (np.std(p)**2) / (np.mean(p)**2) )

    tfinal = max(data[:,3])

pxMF = np.mean( pxM )
pyMF = np.mean( pyM )
pzMF = np.mean( pzM )
pxVF = np.mean( pxV )
pyVF = np.mean( pyV )
pzVF = np.mean( pzV )
pMF  = np.mean( pM )
pVF  = np.mean( pV )
eM   = np.mean( e )

pMS = np.std(pM)
pVS = np.std(pV)
eS  = np.std(e)

with open( 'ptot.dat', 'a') as fp:
    # write out mean polarization and its std.dev and write polarization variance and its std.dev
    s = str(pMF) +'  '+ str(pMS) +'  '+ str(pVF) +'  '+ str(pVS) +'  '+ str(eM) +'  '+ str(eS) +'  '+ str(tfinal) + '\n'
    fp.write(s)

with open( 'pcmp.dat', 'a') as fp:
    s = str(pxMF) +'  '+ str(pyMF) +'  '+ str(pzMF) +'  '+ str(pxVF) +'  '+ str(pyVF) +'  '+ str(pzVF) +'  '+ str(tfinal) + '\n'
    fp.write(s)
