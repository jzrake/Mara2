mesh = GetActiveSource()

slice1 = Slice(mesh)
slice2 = Slice(mesh)

for slice, z in zip([slice1, slice2], [-0.49, 0.49]):
	slice.SliceType = 'Plane'
	slice.SliceType.Origin = [0, 0, z]
	slice.SliceType.Normal = [0, 0, 1]
	slice.Triangulatetheslice = False

	sDisplay = Show(slice)
	sDisplay.Representation = 'Outline'
	sDisplay.ColorArrayName = ['CELLS', '']

	glyph = Glyph(slice)
	glyph.Scalars = ['POINTS', 'None']
	glyph.Vectors = ['CELLS', 'velocity']
	glyph.ScaleMode = 'vector'
	glyph.ScaleFactor = 0.5

	gDisplay = Show(glyph)
	gDisplay.Representation = 'Surface'
	gDisplay.ColorArrayName = ['POINTS', 'velocity']

streamTracer = StreamTracer(mesh)
streamTracer.SeedType.Point1 = [0.0, -0.5, -0.5]
streamTracer.SeedType.Point2 = [0.0, +0.5, -0.5]
streamTracer.MaximumStreamlineLength = 1.0
streamTracer.SeedType.Resolution = 100
Show(streamTracer)
