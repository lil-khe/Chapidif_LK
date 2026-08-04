class NullApp:
    class _NullRunPlot:
        x_axis_keV = True
        def initcalc(self):
            pass
 
    def __init__(self):
        self.runplot = self._NullRunPlot()
 
    def UpdateStatus(self, text=""):
        if text:
            print(text)
 
    def update_all_tables(self):
        pass
