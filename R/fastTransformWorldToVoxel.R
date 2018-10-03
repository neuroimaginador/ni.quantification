## TransformWorldToVoxel 2
fastTransformWorldToVoxel <- function(points, voxelDims) {
  
  new_points <- points
  
  # apply(points, 1, function(x) x/abs(voxelDims)) + 1
  for(i in 1:3) {
    
    new_points[, i] <- new_points[, i]/voxelDims[i]
    
  }
  
  return (new_points + 1)
  
} 