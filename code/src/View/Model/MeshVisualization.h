#pragma once

enum MeshVisualization
{
	Standard,
	LabelAssignment,
	VertexLabelling,
	FaceLabelling,
	PrincipalCurvature,
	GaussianCurvature,
	GaussianCurvaturePrincipal,
	MeanCurvature,
	ErrorDistances
};

struct OverlaySettings
{
	bool has_changed = false;

	float visualized_min = -0.25;
	float visualized_max = 0.25;

	bool show_filtered_curvature_points = false;
	bool show_filtered_curvature_values = false;
	float filter_curvature = 1.0;
};