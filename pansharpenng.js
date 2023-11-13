
/////////////LAND USE AND LAND COVER CLASSIFICATION\\\\\\\
var boundary = Cath_Peak

var s2 = ee.ImageCollection("COPERNICUS/S2_SR")
var rgbVis = {
  min: 0.0,
  max: 3000,
  bands: ['B4', 'B3', 'B2'],
};


// Function to remove cloud and snow pixels from Sentinel-2 SR image

function maskCloudAndShadowsSR(image) {
  var cloudProb = image.select('MSK_CLDPRB');
  var snowProb = image.select('MSK_SNWPRB');
  var cloud = cloudProb.lt(10);
  var scl = image.select('SCL'); 
  var shadow = scl.eq(3); // 3 = cloud shadow
  var cirrus = scl.eq(10); // 10 = cirrus
  // Cloud probability less than 10% or cloud shadow classification
  var mask = cloud.and(cirrus.neq(1)).and(shadow.neq(1));
  return image.updateMask(mask);
}






var image = s2
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5))
  .filter(ee.Filter.date('2019-12-01', '2022-04-30'))
  .filter(ee.Filter.bounds(boundary))
  .map(maskCloudAndShadowsSR)
  //.median()
  .select('B.*')
  
  
var sharpened = image.map(function(image) {
  return panSharpen({
    image: image,
    bestEffort: true
  })
})

print('sharpened', sharpened)
var visParams = {bands: 'B5,B8A,B12', min: 600, max: 4000}
function panSharpen(params) {
  if (params && !(params.image instanceof ee.Image))
    throw Error('panSharpen(params): You must provide an params object with an image key.')
    
  var image = params.image
  var geometry = params.geometry || image.geometry()
  var crs = params.crs || image.select(0).projection()
  var maxPixels = params.maxPixels
  var bestEffort = params.bestEffort || false
  var tileScale = params.tileScale || 1
  
  image = image.clip(geometry)
  
  var bands10m = ['B2', 'B3', 'B4', 'B8']
  var bands20m = ['B5', 'B6', 'B7', 'B8A', 'B11', 'B12']
  var panchromatic = image
    .select(bands10m)
    .reduce(ee.Reducer.mean())
  var image20m = image.select(bands20m)
  var image20mResampled = image20m.resample('bilinear')
  
  var stats20m = image20m
    .reduceRegion({
      reducer: ee.Reducer.stdDev().combine(
        ee.Reducer.mean(), null, true
      ),
      geometry: geometry,
      scale: 20,
      crs: crs, 
      bestEffort: bestEffort, 
      maxPixels: maxPixels, 
      tileScale: tileScale
    })
    .toImage()
    
  var mean20m = stats20m
    .select('.*_mean')
    .regexpRename('(.*)_mean', '$1')
    
  var stdDev20m = stats20m
    .select('.*_stdDev')
    .regexpRename('(.*)_stdDev', '$1')
    
  var kernel = ee.Kernel.fixed({
    width: 5,
    height: 5, 
    weights: [
      [-1, -1, -1, -1, -1],
      [-1, -1, -1, -1, -1],
      [-1, -1, 24, -1, -1],
      [-1, -1, -1, -1, -1],
      [-1, -1, -1, -1, -1]
    ], 
    x: -3, 
    y: -3, 
    normalize: false
  })
  
  var highPassFilter = panchromatic
    .convolve(kernel)
    .rename('highPassFilter')

  var stdDevHighPassFilter = highPassFilter
    .reduceRegion({
      reducer: ee.Reducer.stdDev(),
      geometry: geometry,
      scale: 10,
      crs: crs, 
      bestEffort: bestEffort, 
      maxPixels: maxPixels, 
      tileScale: tileScale
    })
    .getNumber('highPassFilter')
  
  function calculateOutput(bandName) {
    bandName = ee.String(bandName)
    var W = ee.Image().expression(
      'stdDev20m / stdDevHighPassFilter * modulatingFactor', {
        stdDev20m: stdDev20m.select(bandName),
        stdDevHighPassFilter: stdDevHighPassFilter,
        modulatingFactor: 0.25
      }
    )
    return ee.Image()
      .expression(
        'image20mResampled + (HPF * W)', {
          image20mResampled: image20mResampled.select(bandName),
          HPF: highPassFilter,
          W: W
      }
    )
    .uint16()
  }
  
  var output = ee.ImageCollection(
      bands20m.map(calculateOutput)
    )
    .toBands()
    .regexpRename('.*_(.*)', '$1')
    
  var statsOutput = output
    .reduceRegion({
      reducer: ee.Reducer.stdDev().combine(
        ee.Reducer.mean(), null, true
      ),
      geometry: geometry,
      scale: 10,
      crs: crs, 
      bestEffort: bestEffort, 
      maxPixels: maxPixels, 
      tileScale: tileScale
    })
    .toImage()
    
  var meanOutput = statsOutput
    .select('.*_mean')
    .regexpRename('(.*)_mean', '$1')
    
  var stdDevOutput = statsOutput
    .select('.*_stdDev')
    .regexpRename('(.*)_stdDev', '$1')
  
  var sharpened = ee.Image()
    .expression(
      '(output - meanOutput) / stdDevOutput * stdDev20m + mean20m', {
        output: output,
        meanOutput: meanOutput,
        stdDevOutput: stdDevOutput,
        stdDev20m: stdDev20m,
        mean20m: mean20m
      }
    )
    .uint16() 
  
  return image
    .addBands(sharpened, null, true)
    .select(image.bandNames())
}  