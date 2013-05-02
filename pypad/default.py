
# THIS FILE IS PART OF PyPad, AND IS GOVERENED BY A PERMISSIBILITY LICENSE 
# GOVERNING ITS USE AND DISTRIBUTION. YOU SHOULD HAVE RECIEVED A COPY OF THIS
# LICENSE WITH THE SOFTWARE; IF NOT PROVIDED, WRITE TO <tjlane@stanford.edu>.
#
# AUTHORS:
# TJ Lane <tjlane@stanford.edu>
# Jonas Sellberg <jonas.a.sellberg@gmail.com>
#
# Apr 30, 2013

"""
This is a default metrology for the CSPad. It's actually a real metrology:
DSD (DS1) for run 7 at CXI, taken Apr 2013.
"""

import sys
import os
import numpy as np


q0 = np.array( [[0.0, 0.0, 0.0],
                [22.0, 20908.0, 21.0],
                [43523.0, 20942.0, -43.0],
                [43536.0, 40.0, -48.0],
                [12.0, 23318.0, 27.0],
                [-21.0, 44227.0, 26.0],
                [43519.0, 44310.0, -46.0],
                [43556.0, 23397.0, -25.0],
                [271.0, 47051.0, 3.0],
                [278.0, 90590.0, 27.0],
                [21197.0, 90598.0, 36.0],
                [21174.0, 47040.0, 30.0],
                [23756.0, 46994.0, 5.0],
                [23506.0, 90534.0, 31.0],
                [44417.0, 90648.0, 14.0],
                [44661.0, 47088.0, 4.0],
                [47237.0, 45933.0, 3.0],
                [47339.0, 66857.0, -7.0],
                [90876.0, 66638.0, -62.0],
                [90775.0, 45730.0, -63.0],
                [47366.0, 69262.0, 27.0],
                [47425.0, 90172.0, -8.0],
                [90963.0, 90044.0, -30.0],
                [90902.0, 69144.0, -9.0],
                [44494.0, -282.0, -42.0],
                [44291.0, 43266.0, -27.0],
                [65201.0, 43364.0, -27.0],
                [65404.0, -183.0, -52.0],
                [67790.0, -227.0, -60.0],
                [67854.0, 43320.0, -6.0],
                [88760.0, 43289.0, -13.0],
                [88697.0, -264.0, -71.0]] )


q1 = np.array( [[0.0, 0.0, 0.0],
                [14.0, 20919.0, -36.0],
                [43567.0, 20868.0, 30.0],
                [43541.0, -39.0, 70.0],
                [28.0, 23364.0, -26.0],
                [73.0, 44278.0, -87.0],
                [43617.0, 44183.0, -24.0],
                [43565.0, 23265.0, 46.0],
                [350.0, 46662.0, -71.0],
                [444.0, 90582.0, -144.0],
                [21824.0, 90547.0, -90.0],
                [21735.0, 46629.0, -10.0],
                [23737.0, 46548.0, -17.0],
                [23833.0, 90468.0, -87.0],
                [45207.0, 90431.0, -54.0],
                [45118.0, 46504.0, -20.0],
                [47347.0, 45822.0, -11.0],
                [47502.0, 66747.0, 1.0],
                [91045.0, 66420.0, 104.0],
                [90887.0, 45502.0, 79.0],
                [47504.0, 69048.0, -14.0],
                [47557.0, 89963.0, -42.0],
                [91099.0, 89846.0, 47.0],
                [91049.0, 68939.0, 100.0],
                [44377.0, -268.0, 69.0],
                [44413.0, 43274.0, -1.0],
                [65325.0, 43258.0, 53.0],
                [65289.0, -296.0, 134.0],
                [67726.0, -300.0, 143.0],
                [67805.0, 43241.0, 76.0],
                [88713.0, 43213.0, 97.0],
                [88637.0, -338.0, 191.0]] )


q2 = np.array( [[0.0, 0.0, 0.0],
                [-29.0, 20909.0, 8.0],
                [43514.0, 20978.0, -35.0],
                [43533.0, 62.0, 29.0],
                [6.0, 23383.0, 17.0],
                [-11.0, 44298.0, 22.0],
                [43298.0, 44099.0, -2.0],
                [43550.0, 23414.0, 32.0],
                [116.0, 46907.0, 28.0],
                [56.0, 90450.0, 6.0],
                [20955.0, 90478.0, -4.0],
                [21011.0, 46931.0, 30.0],
                [23606.0, 46788.0, 12.0],
                [23389.0, 90329.0, 14.0],
                [44300.0, 90441.0, -17.0],
                [44520.0, 46889.0, 6.0],
                [46962.0, 46079.0, 1.0],
                [46960.0, 67004.0, -6.0],
                [90495.0, 67011.0, -25.0],
                [90506.0, 46098.0, -20.0],
                [47015.0, 69401.0, 9.0],
                [46986.0, 90305.0, -21.0],
                [90526.0, 90357.0, -51.0],
                [90552.0, 69451.0, -25.0],
                [44299.0, 60.0, 15.0],
                [44251.0, 43598.0, 2.0],
                [65159.0, 43619.0, -26.0],
                [65204.0, 79.0, 24.0],
                [67588.0, 41.0, 7.0],
                [67738.0, 43588.0, -11.0],
                [88647.0, 43526.0, -22.0],
                [88509.0, -22.0, 3.0]] )


q3 = np.array( [[0.0, 0.0, 0.0],
                [-53.0, 20908.0, -1.0],
                [43497.0, 21007.0, 14.0],
                [43542.0, 102.0, -18.0],
                [-5.0, 23309.0, 11.0],
                [-84.0, 44216.0, 2.0],
                [43462.0, 44397.0, -35.0],
                [43538.0, 23495.0, -6.0],
                [18.0, 47069.0, -13.0],
                [-241.0, 90615.0, -46.0],
                [20562.0, 90742.0, -19.0],
                [20924.0, 47204.0, -8.0],
                [23419.0, 47173.0, -28.0],
                [23013.0, 90723.0, -76.0],
                [43920.0, 90913.0, -63.0],
                [44336.0, 47370.0, -32.0],
                [47040.0, 45967.0, -23.0],
                [46892.0, 66877.0, -45.0],
                [90431.0, 67173.0, -29.0],
                [90578.0, 46270.0, -65.0],
                [46768.0, 69166.0, -43.0],
                [46611.0, 90151.0, -72.0],
                [90315.0, 90486.0, -65.0],
                [90481.0, 69503.0, -35.0],
                [44377.0, 27.0, -36.0],
                [44231.0, 43574.0, -33.0],
                [65134.0, 43642.0, -31.0],
                [65284.0, 94.0, -26.0],
                [67718.0, 83.0, -32.0],
                [67627.0, 43635.0, -36.0],
                [88536.0, 43669.0, -47.0],
                [88617.0, 123.0, -43.0]] )


