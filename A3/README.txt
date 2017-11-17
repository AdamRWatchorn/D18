As a note before the features, cylinders were attempted but never perfected. This has been a recurring problem, though it is unclear why it has proven to be so. Due to this cylinders are not used and do not have advanced features specifically for them implemented.

Antialiasing was implemented.

Texture Mapping for planes and spheres was implemented. This feature can be seen in the final render.

Area Light Sources for planes was implemented. This can be seen in all renders. Unfortunately spherical area light sources didn't make it into the final submission.

Refractive objects were implemented as can be seen in refraction.ppm. The current implementation of refraction does not support nested refraction. Unfortunately refractive objects didn't make it into the final render, thus refraction.ppm was required, because an unexpected glitch appeared with the objects whenever normal maps are present in the scene. The glitch causes a segfault which failed to be understood and remedied before final submission. The segfault glitch is not present when there is no normal mapping in the scene. Since normal maps were a bonus feature I am hoping that should not be a major issue.

Multi-threading was not implemented.

Hierarchical objects were created in their simplest forms. This is how the brick walls were made in Final_Render.ppm. This was both to have the basics of hierarchical objects in the scene as well as to make moving around the walls easier during development. Unfortunately hierarchical objects of greater complexity and of random nature did not make it into any scene as time grew too short.

Normal Mapping was implemented for planes, which can be seen on every plane in Final_Render.ppm. As mentioned above, normal mapping causes refraction to stop working when it is present. Without normal mapped objects on the scene refractive objects continue to work. Unfortunately spherical normal maps didn't make the final submission.

Alpha mapping, depth-of-field, photon mapping, dispersion, and scene octrees were not implemented for this assignment.
