import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import pandas as pd
import shapely
import matplotlib

print(matplotlib.get_backend())

fig = plt.figure(figsize=(20, 10), frameon=False)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
# ax.axis("off")
ax.set_global()
ax.coastlines()
ax.stock_img()

grid_gdf = gpd.read_file("grid_30deg.gpkg")
grid_gdf = grid_gdf.to_crs(ccrs.Robinson())
# world_gdf = gpd.read_file("world_simplified.gpkg")
# world_gdf = world_gdf.to_crs(ccrs.Robinson())
paralels_gdf = gpd.read_file("parallels_6min.gpkg")
paralels_gdf = paralels_gdf.to_crs(ccrs.Robinson())

grid_gdf.plot(ax=ax, edgecolor="black", linewidth=2.5, zorder=3)
# world_gdf.boundary.plot(ax=ax, linewidth=1.7, color="black", zorder=2)
paralels_gdf.plot(ax=ax, linewidth=0.45, column="ABS Diff (auth inv OLD - common) E-12", cmap="RdYlGn_r", zorder=1)

# ax.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
# ax.plot(grid_gdf.geometry.xy)
# ax.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic())
# plt.subplots_adjust(right=10, left=9, bottom=9, top=10)

plt.savefig("fig.png", transparent=True)

plt.show()
