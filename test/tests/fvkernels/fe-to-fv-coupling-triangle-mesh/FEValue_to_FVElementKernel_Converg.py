#!/usr/bin/env python3

#Uniform element size
import mms
df1 = mms.run_spatial('2D_Coupling_FEValue_to_FVElementKernel.i', 6, 'Mesh/gen/file=triangle_structured_mesh-uniform_progression.msh', y_pp=['u_l2Error', 'v_l2Error'])

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label=['u_l2Error', 'v_l2Error'], marker='o', markersize=8, num_fitted_points=5)
fig.save('FEValue-to-FVElementKernel-uniform-mesh.png')

#Non-uniform element size
import mms
df1 = mms.run_spatial('2D_Coupling_FEValue_to_FVElementKernel.i', 6,'Mesh/gen/file=triangle_structured_mesh-nonuniform_progression.msh', y_pp=['u_l2Error', 'v_l2Error'])

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label=['u_l2Error', 'v_l2Error'], marker='o', markersize=8, num_fitted_points=5)
fig.save('FEValue-to-FVElementKernel-nonuniform-mesh.png')
