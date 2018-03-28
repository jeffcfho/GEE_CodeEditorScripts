//Conversion of Raw DN to spectral radiance + DO subtraction
var convertToLplusDOsubtract = function(image) {
  var Llambda = ee.Algorithms.Landsat.calibratedRadiance(image);
  var geom = image.geometry()
  var DO = Llambda.reduceRegion({
          reducer: ee.Reducer.min(),
          scale:30,
          geometry: geom,
          bestEffort: true});
  return Llambda.subtract(DO.toImage())
    .copyProperties(image,['system:time_end']);
};


var col = ee.ImageCollection('LANDSAT/LT5_L1T').filterBounds(ee.Feature.Point(-83,41.76))
  .filterDate(new Date('7/1/2002'),new Date('10/31/2002'));
var LminusDO = col.map(convertToLplusDOsubtract);
print(col.first());
print(LminusDO.first());

Map.addLayer(ee.Image(col.first()),{bands:['B3', 'B2', 'B1']}, 'Raw DN');
Map.addLayer(ee.Image(LminusDO.first()),{'min': 10, 'max': 60, bands:['B3', 'B2', 'B1']},
  'DO-subtracted spectral radiance ')
Map.setCenter(-82.5,41.8,9);