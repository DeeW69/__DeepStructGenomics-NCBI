import vtkmodules.vtkInteractionStyle  # noqa: F401
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401  (force backend OpenGL2)

from vtkmodules.vtkFiltersSources import vtkCubeSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor


print("VTK: starting...")

src = vtkCubeSource()
src.Update()

mapper = vtkPolyDataMapper()
mapper.SetInputConnection(src.GetOutputPort())

actor = vtkActor()
actor.SetMapper(mapper)

renderer = vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(0.1, 0.1, 0.15)

win = vtkRenderWindow()
win.SetWindowName("VTK Debug Window")
win.SetSize(900, 700)
win.AddRenderer(renderer)

iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(win)

renderer.ResetCamera()
win.Render()

print("VTK: entering interactor.Start() (should block and show a window)")
iren.Initialize()
iren.Start()

print("VTK: interactor ended (window closed)")
