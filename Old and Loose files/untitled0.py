#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""


from __future__ import division
import sympy
m_A, rho_A, rho_B, V, C_pa, C_pb, T= sympy.var('m_A, rho_A, rho_B, V, C_pa, C_pb, T')

func = (1/(m_A * (1+ rho_B/rho_A) - V * rho_B) ) \
        * (m_A * (rho_A + rho_B**2 / rho_A) - V * rho_B**2 ) \
        * (m_A * (C_pa + rho_B*C_pb/ rho_A) - V * rho_B* C_pb ) \
        * T