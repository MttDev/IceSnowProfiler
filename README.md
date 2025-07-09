# SnowIceDepth

**Algorithm for estimating snow and ice thickness using thermistor chain temperature profiles.**

## Overview
This repository contains a simple MATLAB implementation that detects snow and ice layer boundaries by fitting two linear segments to temperature-depth data.

## Contents
- `data.mat` — Example dataset with temperature and depth vectors.
- `profiler.m` — Main script that loads the data, cleans it, applies control charts, and computes the best fitting lines.

## Requirements
- MATLAB (tested with R202x)
- No special toolboxes required.

## Notes
The algorithm was developed for processing thermistor chain data in Arctic conditions.

## License
MIT License.
