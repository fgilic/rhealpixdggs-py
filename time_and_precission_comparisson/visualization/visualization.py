import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import pandas as pd
import shapely
import matplotlib as mpl
from matplotlib.ticker import LogFormatter
import mapclassify
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = plt.figure(figsize=(20, 10), frameon=False)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
# ax.axis("off")
ax.set_global()
ax.coastlines(linewidth=1.75)
# ax.stock_img()

grid_gdf = gpd.read_file("grid_30deg.geojson")
grid_gdf = grid_gdf.to_crs(ccrs.Robinson())
# world_gdf = gpd.read_file("coastlines_simplified.geojson")
# world_gdf = world_gdf.to_crs(ccrs.Robinson())
paralels_gdf = gpd.read_file("parallels_6min.geojson")
paralels_gdf = paralels_gdf.to_crs(ccrs.Robinson())

grid_gdf.plot(ax=ax, edgecolor="black", linewidth=2.5, zorder=3)
# world_gdf.boundary.plot(ax=ax, linewidth=1.1, color="black", zorder=2)

# https://gis.stackexchange.com/questions/330008/center-normalize-choropleth-colors-in-geopandas
# norm = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=0.6, vmax=5)

# # abs_diff_auth_old-mpmath.svg
# norm = mpl.colors.PowerNorm(gamma=0.12, vmin=0, vmax=5)
# paralels_gdf.iloc[:-1].plot(ax=ax, linewidth=0.35, column="ABS Diff (auth OLD - auth direct mpmath) E-12",
#                             cmap="RdYlGn_r", zorder=1, legend=True, norm=norm)


# geodetic -> authalic
old_new_concat = pd.concat(
    [
        paralels_gdf["ABS Diff (auth NEW - auth direct mpmath) E-12"],
        paralels_gdf["ABS Diff (auth OLD - auth direct mpmath) E-12"],
    ]
).iloc[:-1]
vmin = old_new_concat.min()
vmax = old_new_concat.max()
# values in (value, color) defined by bin edges with Quantiles, 6 bins
# mapclassify.classify(old_new_concat, scheme="Quantiles", k=6)
# 0, 0.0013390670000000004, 0.003919687, 0.007541177, 0.014841380000000038, 0.03943446300000004, 8.869908902
# and normalized from 0.0 to 1.0
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    N=50000000,
    name="geodetic-auth",
    colors=[
        (0, "#006837"),
        (0.000150967, "#4bb05c"),
        (0.000441908, "#b7e075"),
        (0.000850198, "#fffebe"),
        (0.001673228, "#fdbf6f"),
        (0.00444587, "#ea5739"),
        (1, "#a50026"),
    ],
)

# # plot errors obtained by current statements (geodetic -> authalic)
# paralels_gdf.iloc[:-1].plot(
#     ax=ax,
#     linewidth=0.35,
#     column="ABS Diff (auth OLD - auth direct mpmath) E-12",
#     cmap=cmap,
#     zorder=1,
#     legend=True,
#     vmin=vmin,
#     vmax=vmax,
# )
# plt.savefig("geodetic-authalic_current.svg")

# # plot errors obtained by modified statements (geodetic -> authalic)
# paralels_gdf.iloc[:-1].plot(
#     ax=ax,
#     linewidth=0.35,
#     column="ABS Diff (auth NEW - auth direct mpmath) E-12",
#     cmap=cmap,
#     zorder=1,
#     legend=True,
#     vmin=vmin,
#     vmax=vmax,
# )
# plt.savefig("geodetic-authalic_modified.svg")


# authalic -> geodetic
old_new_concat = pd.concat(
    [
        paralels_gdf["ABS Diff (auth inv OLD - common) E-12"],
        paralels_gdf["ABS Diff (auth inv NEW - common) E-12"],
    ]
)
vmin = old_new_concat.min()
vmax = old_new_concat.max()
# values in (value, color) defined by bin edges with Quantiles, 6 bins
# mapclassify.classify(old_new_concat, scheme="Quantiles", k=6)
# 0, 0, 0.004000, 0.030000, 1405.600000, 7976.770000, 14201.336000
# and normalized from 0.0 to 1.0
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    N=50000000,
    name="auth-geodetic",
    colors=[
        (0, "#006837"),
        (0, "#4bb05c"),
        (2.82e-07, "#b7e075"),
        (2.11e-06, "#fffebe"),
        (0.098977, "#fdbf6f"),
        (0.561692, "#ea5739"),
        (1, "#a50026"),
    ],
)
# # plot errors obtained by current statements (authalic -> geodetic)
# paralels_gdf.plot(
#     ax=ax,
#     linewidth=0.35,
#     column="ABS Diff (auth inv OLD - common) E-12",
#     cmap=cmap,
#     zorder=1,
#     legend=True,
#     vmin=vmin,
#     vmax=vmax,
# )
# plt.savefig("authalic-geodetic_current.svg")

# plot errors obtained by current statements (authalic -> geodetic)
paralels_gdf.plot(
    ax=ax,
    linewidth=0.35,
    column="ABS Diff (auth inv NEW - common) E-12",
    cmap=cmap,
    zorder=1,
    legend=True,
    vmin=vmin,
    vmax=vmax,
)
plt.savefig("authalic-geodetic_modified.svg")


# ax.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
# ax.plot(grid_gdf.geometry.xy)
# ax.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic())
# plt.subplots_adjust(right=10, left=9, bottom=9, top=10)

plt.show()
