class NullApp:
    """
    Stand-in for the Qt main window that `calculate` (and, transitively,
    `run_and_plot`) expects as `parent`. Covers every self.MyChapApp.*
    call site currently in calculate.py:
        - UpdateStatus(text)      -> just prints, so you get progress
                                      messages in the console instead of
                                      a GUI status bar
        - update_all_tables()     -> no-op
        - runplot.x_axis_keV      -> plain attribute, defaults to True
        - runplot.initcalc()      -> no-op
    If you hit a fresh AttributeError for some other self.MyChapApp.X,
    it just means a call site I haven't seen yet -- add X here the same way.
    """
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
