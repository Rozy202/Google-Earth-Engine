// Create a polygon feature using the provided coordinates for "IX"
var iiiPolygon = ee.Geometry.Polygon([
  [
    [29.236167408325155, -29.00419162733362],
    [29.236194188262104, -29.003536103595934],
    [29.23622087951898, -29.002671008273065],
    [29.236483978795448, -29.001859492749716],
    [29.23674264386018, -29.001128201503434],
    [29.236978962425006, -29.000731293893587],
    [29.236978962425006, -29.00007580216075],
    [29.23750515783917, -28.99926425425955],
    [29.237826159185424, -28.998234252850878],
    [29.23800900549973, -28.9979666908209],
    [29.238062564707434, -28.99762780927361],
    [29.238231932274076, -28.997181862307293],
    [29.238396954623877, -28.996851903417813],
    [29.238633270567124, -28.996120644930222],
    [29.238762556954683, -28.995411629235857],
    [29.238816115891716, -28.994702611581875],
    [29.239079210589484, -28.9939177949711],
    [29.23957418560447, -28.993025970912623],
    [29.239966564308013, -28.992241232374756],
    [29.240412500342366, -28.991612415920635],
    [29.24114821651111, -28.991113007941006],
    [29.241224031337183, -28.991072913298407],
    [29.241415740701953, -28.99077418053498],
    [29.24092529409737, -28.98962812629797],
    [29.240269737112047, -28.988928054193273],
    [29.23965443489225, -28.988437572275238],
    [29.238753689579134, -28.987987207554742],
    [29.23759879808228, -28.987808843170104],
    [29.23708155914877, -28.98798719312481],
    [29.23678272493395, -28.988187879430722],
    [29.235926654702574, -28.98852231949012],
    [29.234999106726764, -28.989151005893014],
    [29.234316832096447, -28.989409712357343],
    [29.23395565269856, -28.98954341461752],
    [29.233179819015078, -28.990323766122625],
    [29.232279034601596, -28.991674907099863],
    [29.23174395296148, -28.99237052809253],
    [29.23129364577904, -28.993355955391294],
    [29.2308833312181, -28.99454206695736],
    [29.230473104207874, -28.995242147082614],
    [29.230107481562197, -28.996263339312886],
    [29.229572392137026, -28.99802908838656],
    [29.228957045794942, -28.999745864194516],
    [29.22859141943614, -29.00114156184981],
    [29.22797606905088, -29.002287565717655],
    [29.226745360672584, -29.00343796967652],
    [29.225630548990722, -29.004441294361282],
    [29.226879091401166, -29.004441296639786],
    [29.228158838755615, -29.004311954887907],
    [29.230334854081637, -29.005145813340786],
    [29.23282750459318, -29.005145818184847],
    [29.234428386960854, -29.004887228820117],
    [29.23602472970586, -29.00437889929335],
    [29.236207578228775, -29.004365452512552],
    [29.236167408325155, -29.00419162733362]
  ]
]);



var SMAP = ee.ImageCollection('NASA_USDA/HSL/SMAP10KM_soil_moisture')
  .filterBounds(Cath_Peak)
  .filterDate('2018-11-14', '2022-07-29');

var proj_10m = ee.Image.pixelLonLat().projection();

var soilmoisture = SMAP.select('smp');

function resample(image) {
  var soilmoisture = image.select('smp');
  var soilmoisture_res2 = soilmoisture.resample('bilinear').reproject(proj_10m);
  return soilmoisture_res2;
}

var soilmoisture_resampled = soilmoisture.map(resample);

// Function to extract data for the "IX" polygon within Cath Peak
var extractData = function(image) {
  var value = image.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: ixPolygon,
    scale: 100 // Adjust the scale to control the output resolution
  });
  return ee.Feature(null, value)
    .set('system:time_start', image.get('system:time_start'));
};

// Map over the ImageCollection to extract data for the "IX" polygon
var data = soilmoisture_resampled.map(extractData);

// Export the data to Google Drive as a CSV file
Export.table.toDrive({
  collection: ee.FeatureCollection(data),
  description: 'SMP_catch3_Export',
  fileFormat: 'CSV'
});