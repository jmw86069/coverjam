
# coverjam 0.0.1.900

* Migrated the updated `nmatlist2heatmaps()` from platjam.
* Migrated `zoom_nmat()` and `zoom_nmatlist()` from platjam.
* `nmatlist2heatmaps()`

   * default `raster_by_magick` checks whether magick is installed.
   While `magick` is hugely preferred, it is not mandatory for `coverjam`
   because it can be non-trivial to install on some linux systems,
   mainly by requiring development libraries which may require an updated
   build chain, which ultimately requires root.
   * added `magick` to Suggests.

# coverjam 0.0.0.900

Initial commit. Ported the first set of functions from `platjam`.

