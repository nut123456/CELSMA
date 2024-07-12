# -*- coding: utf-8 -*-
"""
Created on Fri May 31 14:42:09 2024

@author: PC
"""

import matlab.engine
import numpy as np

# เริ่มต้น MATLAB Engine
eng = matlab.engine.start_matlab()

# ระบุโฟลเดอร์ที่มีฟังก์ชัน MATLAB
eng.addpath(r'C:\Users\PC\OneDrive - Chulalongkorn University\Desktop\PAPER\Dome\600_Bar_benchmark(matlab)')

# กำหนดค่าอินพุต x
x = np.random.rand(25).tolist()

# เรียกใช้ฟังก์ชัน Bar600_truss ใน MATLAB
W, f2 = eng.Bar600_truss(x, nargout=2)
print(f"W: {W}")
print(f"f2: {f2}")

# ปิด MATLAB Engine
eng.quit()
