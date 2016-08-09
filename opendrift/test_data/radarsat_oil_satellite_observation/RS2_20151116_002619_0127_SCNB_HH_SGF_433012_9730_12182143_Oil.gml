<?xml version='1.0' encoding='ISO-8859-1'?>
<od:oilDetectionCollection xmlns:od="http://cweb.ksat.no/cweb/schema/geoweb/oil" xmlns:gml="http://www.opengis.net/gml" xmlns:xsi="http://www.w4.org/2001/XMLSchema-instance" xsi:schemaLocation="http://cweb.ksat.no/cweb/schema/geoweb/oil http://cweb.ksat.no/cweb/schema/geoweb/oil/oilDetection.xsd">
  <od:oilDetectionMember>
    <od:oilDetection>
      <od:oilSpill>
        <gml:Polygon>
          <gml:exterior>
            <gml:LinearRing>
              <gml:posList srsName="CRS:84">4.043959 60.539858 4.047977 60.540271 4.047561 60.537130 4.055878 60.533165 4.061807 60.531537 4.065977 60.530761 4.075421 60.531558 4.086008 60.532128 4.093012 60.534567 4.095190 60.536855 4.100212 60.537370 4.101578 60.535789 4.100083 60.532710 4.102131 60.530338 4.099605 60.528875 4.101262 60.525774 4.104075 60.521588 4.103054 60.518213 4.098939 60.515554 4.095910 60.511974 4.093930 60.510050 4.092840 60.506840 4.099306 60.504750 4.108724 60.502273 4.113742 60.502787 4.116948 60.503459 4.116632 60.506697 4.119397 60.510078 4.121430 60.508565 4.126082 60.506631 4.131726 60.505660 4.135683 60.505376 4.136169 60.504221 4.136552 60.500819 4.139069 60.498151 4.140916 60.495414 4.142037 60.491915 4.143423 60.488615 4.139030 60.486617 4.132115 60.483156 4.129143 60.480271 4.125099 60.477448 4.110265 60.475412 4.100825 60.475476 4.095537 60.479753 4.087252 60.482861 4.084412 60.484634 4.088105 60.488283 4.084595 60.489988 4.078229 60.490195 4.075400 60.486979 4.067502 60.486683 4.066819 60.491604 4.061274 60.494821 4.065746 60.495797 4.072183 60.495426 4.075528 60.499900 4.067696 60.503571 4.068044 60.506877 4.066287 60.511859 4.062077 60.515213 4.058604 60.518471 4.055116 60.522588 4.050149 60.526896 4.045780 60.531439 4.043098 60.535294 4.043959 60.539858</gml:posList>
            </gml:LinearRing>
          </gml:exterior>
        </gml:Polygon>
      </od:oilSpill>
      <od:guid>4d2b55f4-3689-49f4-80c6-10b89a566919</od:guid>
      <od:oilSpillID>1</od:oilSpillID>
      <od:alertForm>NOTIFICATION</od:alertForm>
      <od:oilDefinition>SPILL</od:oilDefinition>
      <od:detectionTime>2015-11-16T00:26:18.770Z</od:detectionTime>
      <od:geometricMeasurementsProperty>
        <od:geometricMeasurements>
          <od:area>17478931.5921</od:area>
          <od:length>8124.57982169</od:length>
          <od:width>6035.38642502</od:width>
          <od:center>
            <gml:Point>
              <gml:pos srsName="CRS:84">4.094639 60.504534</gml:pos>
            </gml:Point>
          </od:center>
          <od:orientation>299.681579536</od:orientation>
        </od:geometricMeasurements>
      </od:geometricMeasurementsProperty>
      <od:geometricMeasurementsByPatchesProperty>
        <od:geometricMeasurementsByPatches patchIndex="1">
          <od:area>17478931.5921</od:area>
          <od:length>8124.57982169</od:length>
          <od:width>6035.38642502</od:width>
          <od:center>
            <gml:Point>
              <gml:pos srsName="CRS:84">3.094639 60.504534</gml:pos>
            </gml:Point>
          </od:center>
          <od:orientation>299.681579536</od:orientation>
        </od:geometricMeasurementsByPatches>
      </od:geometricMeasurementsByPatchesProperty>
      <od:confidenceCategory>A</od:confidenceCategory>
      <od:classificationProperty>
        <od:classification>
          <od:shapeOutline absRatio="0.0578035">CONTINIOUS</od:shapeOutline>
          <od:shapeIsPatch absRatio="0.0289017">true</od:shapeIsPatch>
          <od:contrastCharacteristic absRatio="0.0578035">VARIABLE</od:contrastCharacteristic>
          <od:edgeCharacteristic absRatio="0.0867052">DIFFUSE</od:edgeCharacteristic>
          <od:textureCharacteristic absRatio="0.0289017">VARIABLE</od:textureCharacteristic>
          <od:possibleSourceAssessment absRatio="0.289017">LIKELY SOURCE</od:possibleSourceAssessment>
          <od:shapeRelatedToWindHistory absRatio="0.0578035">MATCH</od:shapeRelatedToWindHistory>
          <od:isRepeatedObservation absRatio="0">false</od:isRepeatedObservation>
          <od:isNaturalSlicksInVicinity absRatio="0">false</od:isNaturalSlicksInVicinity>
          <od:windSpeed absRatio="0.0462428">2.82168</od:windSpeed>
        </od:classification>
      </od:classificationProperty>
      <od:metocProperty>
        <od:metocWind>
          <od:metocObservation>remote</od:metocObservation>
          <od:metocSource>SARTool - SAR ocean wind</od:metocSource>
          <od:windSpeed>2.82168316841</od:windSpeed>
          <od:windDirection>187.896606445</od:windDirection>
        </od:metocWind>
        <od:metocWind>
          <od:metocObservation>model</od:metocObservation>
          <od:metocSource>ncep - Assimilated 10m wind</od:metocSource>
          <od:windSpeed>8.51639595134</od:windSpeed>
          <od:windDirection>188.16815999</od:windDirection>
        </od:metocWind>
        <od:metocWave>
          <od:metocObservation>model</od:metocObservation>
          <od:metocSource>met.no - Assimilated ocean wave</od:metocSource>
          <od:waveHeight>2.3</od:waveHeight>
          <od:waveDirection>8.5</od:waveDirection>
        </od:metocWave>
      </od:metocProperty>
      <od:externalDataProperty>
        <od:externalRef>
          <od:ref>RS2_20151102_172619_0127_SCNB_HH_SGF_433012_9730_12182143_OSN_dc06963e-b725-4127-89ad-c458fab57df5_4d2b55f4-3689-49f4-80c6-10b89a566919.png</od:ref>
          <od:refName>OILIMAGE</od:refName>
          <od:refDescription>file path to an image subset of the oil spill</od:refDescription>
        </od:externalRef>
      </od:externalDataProperty>
      <od:description> TA-1, P-MEDIUM, K-MEDIUM</od:description>
      <od:productMetadataProperty>
        <od:productMetadata>
          <od:guid>a0fc6364-ce66-4248-a984-343c90288681</od:guid>
          <od:productID>RS2_20151102_172619_0127_SCNB_HH_SGF_433012_9730_12182143</od:productID>
          <od:startTime>2015-11-02T17:26:18.770Z</od:startTime>
          <od:stopTime>2015-11-02T17:28:25.315Z</od:stopTime>
        </od:productMetadata>
      </od:productMetadataProperty>
      <od:possibleSourceProperty>
        <od:possibleSource>
          <od:objectLocation>
            <gml:Point>
              <gml:pos srsName="CRS:84">3.046006 60.542291</gml:pos>
            </gml:Point>
          </od:objectLocation>
          <od:objectName>BRAGE</od:objectName>
          <od:guid>884d0bde-c716-4de1-8649-b75ecf2dfd2c</od:guid>
          <od:objectType>PLATFORM</od:objectType>
          <od:confidence>1.0</od:confidence>
        </od:possibleSource>
      </od:possibleSourceProperty>
      <od:previousDetectionProperty/>
    </od:oilDetection>
  </od:oilDetectionMember>
</od:oilDetectionCollection>
