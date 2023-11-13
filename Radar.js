////////////////////////////////////////////////////////////////////Radar\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

var collection_S1 = ee.ImageCollection('COPERNICUS/S1_GRD')
.filterDate ('2021-11-01','2022-04-30')
.filterBounds(Cath_Peak)
print('S1 images over Drakensberg', collection_S1);

//FILTER OUT DIFFERENT MODES
var collection_S1_IW = collection_S1.filterMetadata('instrumentMode', 'equals', 'IW')
print('S1 IW mode images',collection_S1_IW)



//ORDER IMAGES IN TIME AND DISPLAY FIRST IMAGE AS LAYER
var collection_S1_IW_first = ee.Image(collection_S1_IW.sort('system:time_start',false).first());


var collection_S1_IW_ASC =collection_S1_IW.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
print('S1 IW images from Ascending orbit',collection_S1_IW_ASC)

var collection_S1_IW_DES = collection_S1_IW.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));


//CREATE A COMPOSITE FROM MEANS AT DIFFERENT POLARIZATIONS AND LOOK ANGLES AND DISPLAY
var collection_S1_IW_composite =ee.Image.cat([
  collection_S1_IW_ASC.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).select('VH').mean(),
  collection_S1_IW_ASC.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).select('VV').mean(),
  collection_S1_IW_DES.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).select('VV').mean(),
  collection_S1_IW_DES.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).select('VH').mean(),
  ]);

// filter to get images from ascending mode and VH
var collection_S1_IW_ASC = collection_S1_IW.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
var collection_S1_IW_ASC_VH = collection_S1_IW_ASC.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).select('VH')

// define function to convert from dB to power and vice versa
function topow(image) {return ee.Image(10.0).pow(image.divide(10.0))}
function todB(image) {return ee.Image(image).log10().multiply(10.0)}


var collection_S1_IW_ASC_VH_std_pow = collection_S1_IW_ASC_VH.map(topow).reduce(ee.Reducer.stdDev())
var collection_S1_IW_ASC_VH_std_dB = collection_S1_IW_ASC_VH.reduce(ee.Reducer.stdDev())
var collection_S1_IW_first_std = collection_S1_IW_first.select('VH').reduceNeighborhood({
  reducer: ee.Reducer.stdDev(),
  kernel: ee.Kernel.circle(5) //5 is the kernel radius in pixels
  })

// apply a Refined Lee filter (as coded in the SNAP 3.0 S1TBX:) and a simple boxcar filter and observe the differences
// The RL speckle filter
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result);
}

// First add the unfiltered image

// Generate a dB version of the RL filtered image stack
var collection_S1_IW_first_VH_RL = todB(RefinedLee(topow(collection_S1_IW_first.select('VH'))))
var collection_S1_IW_first_VV_RL = todB(RefinedLee(topow(collection_S1_IW_first.select('VV'))))

// Compare to simple 3x3 boxcar filtered stack
var collection_S1_IW_first_VH_BC = todB(topow(collection_S1_IW_first.select('VH')).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.rectangle(3,3)))
var collection_S1_IW_first_VV_BC = todB(topow(collection_S1_IW_first.select('VV')).reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.rectangle(3,3)))
//Map.addLayer(collection_S1_IW_first_VV_BC.addBands(collection_S1_IW_first_VH_BC).addBands(collection_S1_IW_first_VV_BC), {min: [-15, -25, -15],max: [0, -5, 0]}, 'Sentinel-1 IW VV,VH,VV Boxcar 3x3 filtered',false);
