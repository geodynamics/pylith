import numpy
import matplotlib.pyplot as pyplot

dof = numpy.array([8*16, 16*32, 32*64, 64*128, 128*256, 256*512, 512*1024])
i2mult = numpy.array([
    [21, 2.5],
    [29, 2.4],
    [38, 2.5],
    [52, 3.2],
    [81, 6.3],
    [119, 19.9],
    [172, 95.3],
    ])
i3lower = numpy.array([
    [30, 2.6],
    [47, 2.7],
    [70, 2.8],
    [98, 5.1],
    [156, 8.3],
    [239, 31.1],
    [391, 160.7],
    ])
i2user = numpy.array([ # rtol=1.0e-4
    [17+7, 3.6],
    [32, 2.6],
    [35, 3.0],
    [36, 5.0],
    [46, 18.5],
    [46, 80.3],
    [49, 445.1],
    ])
i3user = numpy.array([ # rtol=1.0e-2
    [38, 2.9],
    [54, 2.7],
    [72, 3.4],
    [73, 6.1],
    [76, 28.0],
    [79, 132.6],
    [90, 753.1],
])

class Figure():

    def __init__(self):
        self.fig, (self.axL, self.axR) = pyplot.subplots(2, 1)

    def plot(self, dof, miter, rtime, label):
        mdof = dof[:len(miter)]
        ll = self.axL.loglog(mdof, miter, label=label)
        pfit = self.fit(mdof, miter)

        self.axL.loglog(dof, 10**(numpy.log10(dof)*pfit[0] + pfit[1]), '--', color=ll[0].get_color())
        self.axL.set_xlabel("Problem Size (# cells)")
        self.axL.set_ylabel("# KSP Iterations")
    
        self.axR.loglog(mdof, rtime, label=label)
        self.axR.set_xlabel("Problem Size (# cells)")
        self.axR.set_ylabel("Runtime (seconds)")
    
    def fit(self, mdof, miter):
        return numpy.polyfit(numpy.log10(mdof[:len(miter)]), numpy.log10(miter), 1)

    def run(self):
        self.plot(dof, i3user[:,0], i3user[:,1], "v3 Schur upper + user")
        self.plot(dof, i2user[:,0], i2user[:,1], "v2 Schur upper + user")
        self.plot(dof, i2mult[:,0], i2mult[:,1], "v2 Custom PC w/mult")
        self.plot(dof, i3lower[:,0], i3lower[:,1], "v3 Schur lower + selfp")

        self.axR.legend(loc="upper left")
        pyplot.show()

Figure().run()
