from pygimli.physics.SIP import SIPSpectrum
import pandas as pd


# example for fitting by double cole-cole model
loadPhi = pd.read_excel("example for fit.xlsx", sheet_name='phase')
loadCon = pd.read_excel("example for fit.xlsx", sheet_name='conductivity')
fre = pd.DataFrame(loadPhi.iloc[2:39, 0])  # Read x coordinate value: frequency
fre = fre.reset_index(drop=True)  # reset the index from 0
for i in range(0, 3):
    print('fit No.', i)
    phi = pd.DataFrame(loadPhi.iloc[2:39, i + 2],)  # Read the i-th phase value
    phi = phi.reset_index(drop=True) # reset the index from 0
    con = pd.DataFrame(loadCon.iloc[2:39, i + 2], )  # Read the i-th conductivity value
    con = con.reset_index(drop=True)  # reset the index from 0
    sipFit = pd.concat([fre, con, phi], axis=1)  # Get a set of SIP data for fitting
    sipFit.to_csv('sip-variation.csv', header=False, index=False)
    sip = SIPSpectrum('sip-variation.csv', unify=True)
    sip.fitamp()
    sip.showAll()

    mCC, Rho0, phiCC, ampCC = sip.savedata2() # get the corresponding parameters in frequency i
    mCC = pd.DataFrame(mCC)
    Rho0 = pd.DataFrame(Rho0)
    if i == 0:
        par = mCC
        Rho = Rho0
        re = pd.DataFrame(ampCC * np.cos(phiCC))
        im = pd.DataFrame(ampCC * np.sin(phiCC))
    else:
        par = np.append(par, mCC, axis=1) # combine all the mCC (in ndarray)
        Rho = np.append(Rho, Rho0, axis=1)


par = pd.DataFrame(par)
Rho = pd.DataFrame(Rho)  # sigma0 at low frequency
phiCC = pd.DataFrame(phiCC)
with pd.ExcelWriter('fitResults.xlsx') as writer:
    par.to_excel(writer, sheet_name='paramaters', index=False)
    Rho.to_excel(writer, sheet_name='Rho-vary', index=False)
    phiCC.to_excel(writer, sheet_name='phase', index=False)
    im.to_excel(writer, sheet_name='imaginary', index=False)
