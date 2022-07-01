# The double cole-cole model of SIP based on pyGIMLi
This model is modified by open source libray pyGIMLi and used for simulating a SIP signals with double peaks. Here are the steps for guide:
## Preparation
  * Install [pyGIMLi](https://www.pygimli.org/installation.html) according to the official website steps
  * Download [SIPModel](https://github.com/bluspi/The-double-cole-cole-model-of-SIP-based-on-pygimli/tree/main/SIPModel) folder to any location
## Step 1: Add fitamp function in `sipspectrum.py`
```
    def fitamp(self, eRho=0.01, ePhi=0.001, lam=1000., mpar=(0, 0, 1),
                  taupar1=(0, 1e-5, 1), taupar2=(0, 1e-1, 1000),
                  cpar=(0.5, 0, 1), verbose=False):
        '''made by caojingjing'''
        if taupar1[0] == 0:
            taupar1 = (np.sqrt(taupar1[1] * taupar1[2]), taupar1[1], taupar1[2])
            print("taupar1", taupar1)
        if taupar2[0] == 0:
            taupar2 = (np.sqrt(taupar2[1] * taupar2[2]), taupar2[1], taupar2[2])
            print("taupar2", taupar2)
        #            taupar1[0] = 1.0 / self.f[np.argmax(self.phi)] / 2.0 / pi
        if mpar[0] == 0:
            mpar = (1. - min(self.amp) / max(self.amp), mpar[1], mpar[2])  # *2
            print("mpar", mpar)

        model, self.ampCC, self.phiCC = fitdoubleCC(self.f, self.amp, self.phi,
                                         mpar1=mpar, mpar2=mpar,
                                         cpar1=cpar, cpar2=cpar,
                                         taupar1=taupar1, taupar2=taupar2)
        # re, im = self.realimag(cond=cond)
        self.Rho0 = model[:1]
        self.mCC = model[1:]
```
## Step 2: Add fitdoubleCC function in `tools.py`
```
def fitdoubleCCC(f, amp, phi, eRho=0.01, ePhi=0.001, lam=1000.,
                 mpar1=(0.2, 0, 1), taupar1=(1e-2, 1e-5, 100), cpar1=(0.2, 0, 1),
                 mpar2=(0.2, 0, 1), taupar2=(1e-2, 1e-5, 100), cpar2=(0.2, 0, 1)):
    """Fit double C-C amplitude and phase by caojingjing."""
    f2CC = doubleColeCole(f)
    tLog = pg.trans.TransLog()
    f2CC.region(0).setStartModel(max(amp))
    f2CC.region(1).setParameters(*mpar1)    # m (start,lower,upper)
    f2CC.region(2).setParameters(*taupar1)  # tau
    f2CC.region(3).setParameters(*cpar1)  # c
    f2CC.region(4).setParameters(*mpar2)  # m (start,lower,upper)
    f2CC.region(5).setParameters(*taupar2)  # tau
    f2CC.region(6).setParameters(*cpar2)  # c
    data = pg.cat(amp, phi)
    ICC = pg.core.Inversion(data, f2CC, False)  # set up inversion class
    ICC.setTransModel(tLog)
    error = pg.cat(eRho*amp, pg.Vector(len(f), ePhi))
    ICC.setAbsoluteError(error)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = np.asarray(ICC.run())  # run inversion
    ICC.echoStatus()
    response = np.asarray(ICC.response())
    amp, fre = response[:len(f)], response[len(f):]
    return model, amp, fre
   ```
## Step 3: Add doubleColeCole class in `models.py`
```
class doubleColeCole(pg.core.ModellingBase):
    """fit double on complex conductivity by cole-cole model-CJJ"""

    def __init__(self, f, verbose=False):
        super(doubleColeCole, self).__init__(verbose)
        self.f_ = - f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 7))  # 7 single parameters

    def response(self, par):
        """phase angle of the model"""
        spec1 = modelColeColeSigma(self.f_, 1, par[1], par[2], par[3])
        spec2 = modelColeColeSigma(self.f_, 1, par[4], par[5], par[6])
        spec = par[0] * (spec1 + spec2 - 1.)
        a = np.abs(spec)
        b = -np.angle(spec)
        c = pg.cat(np.abs(spec), -np.angle(spec))
        return pg.cat(np.abs(spec), -np.angle(spec))
```
## Step 4: Start your simulation
Run `fitting a double cole-cole model.py`
