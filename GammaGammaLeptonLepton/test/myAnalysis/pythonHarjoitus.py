from ROOT import gROOT, TCanvas, TF1
import time

gROOT.Reset()

c1 = TCanvas('c1','Example',200,10,700,500)

fun = TF1('fun1','abs(sin(x)/x)',0,10)
c1.SetGridx()
c1.SetGridy()
fun.Draw()
c1.Update()

time.sleep(3)
