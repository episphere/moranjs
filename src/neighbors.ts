import { area, polygonContains } from "d3"
import { FeatureCollection, Position, Feature } from "geojson"

import * as d3 from "https://cdn.jsdelivr.net/npm/d3-array@3/+esm"
import * as kd from 'https://cdn.skypack.dev/kd-tree-javascript@1.0.3?min'


export type Method = "rook" | "queen"
export function findNeighbors(featureCollection: FeatureCollection, method:Method="queen") {
  if (method == "rook") {
    return neighborsRook(featureCollection)
  } else {
    const rook = neighborsRook(featureCollection)
    const point = neighborsPoint(featureCollection)
    const map = d3.index(rook, d => d.join("_&_"))
    point.forEach(pair => map.set(pair.join("_&_"), pair))
    return [...map.values()]
  }
}

function featurePoints(feature: Feature, alternativeId: ID) {
  const areaPoints: AreaPoints[] = []

  const id = feature.id != null ? feature.id : alternativeId

  if (feature.geometry.type == "MultiPolygon") {

    feature.geometry.coordinates.forEach(polygon => {
      for (const points of polygon) {
        const bbox = pointsBoundingBox(points)
        if (bbox != null) {
          areaPoints.push({id: id, points: points, bbox: bbox})
        }
      }
    })
  } else if (feature.geometry.type == "Polygon") {
    for (const points of feature.geometry.coordinates) {
      const bbox = pointsBoundingBox(points)
      if (bbox != null) {
        areaPoints.push({id: id, points: points, bbox: bbox})
      }
    }
  }

  return areaPoints
}

// ######## Point Neighbors ###########

function neighborsPoint(featureCollection: FeatureCollection) {
  const allAreaPoints: AreaPoints[] = [] 

  featureCollection.features.forEach((feature, i) => {
    const areaPoints = featurePoints(feature, i)
    areaPoints.forEach(d => allAreaPoints.push(d))
  })

  const overlapping = findOverlapping(allAreaPoints, 0.00001)
  return overlapping
}

type TreePoint = {x: number, y: number, id: ID}
function findOverlapping(allAreaPoints: AreaPoints[], threshold=1) {
  const treePoints:TreePoint[] = []
  allAreaPoints.forEach((areaPoints) => {
    areaPoints.points.forEach(p => treePoints.push({x: p[0], y: p[1], id: areaPoints.id}))
  })

  const tree = new kd.kdTree(treePoints, 
    (a,b) => Math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2), ["x",  "y"])

  const overlapsMap = new Map<ID, Set<ID>>()
  let startSearchN = 2

  treePoints.forEach((point, i) => {
    let existingOverlaps = overlapsMap.get(point.id)
    if (!existingOverlaps) {
      existingOverlaps = new Set()
      overlapsMap.set(point.id, existingOverlaps)
    }
    let overlaps:Set<ID> = existingOverlaps
    
    let searchN = startSearchN
    while (true) {
      const neighbors = tree.nearest(point, searchN, threshold)

    
      if (neighbors.length < searchN){ 
        neighbors.filter(neighbor => neighbor[0].id != point.id).forEach(d => overlaps.add(d[0].id))
        break 
      } else {
        searchN = searchN * 2
      }
    }

    startSearchN = Math.max(1, Math.ceil(startSearchN + (overlaps.size - startSearchN)/(i+1)))
  })
  
  const overlapPairs:[ID,ID][] = [] 
  for (const [k,v] of overlapsMap.entries()) {
    v.forEach(d => overlapPairs.push([k,d]))
  }

  return overlapPairs
}

// ######## Rook Neighbors ###########

export type ID = string | number
type Arc = [Position, Position]
type LineEquation = {m: number, c: number}
type ArcLine = {points: Arc, line: LineEquation}
type AreaPoints = {id: ID, points: Position[], bbox: [Position, Position]}
type AreaArcs = {id: ID, arcs: ArcLine[]}
function neighborsRook(featureCollection: FeatureCollection):[ID,ID][] {
  const areaPoints: AreaPoints[] = [] 

  featureCollection.features.forEach((feature, i) => {
    const points = featurePoints(feature, i)
    points.forEach(d => areaPoints.push(d))
  })

  const possibleNeighborPolygons:[number,number][] = []
  for (let i = 0; i < areaPoints.length-1; i++) {
    for (let j = i+1; j < areaPoints.length; j++) {
      if (rectanglesOverlap(areaPoints[i].bbox, areaPoints[j].bbox, 5)) {
        possibleNeighborPolygons.push([i,j])
      }
    }
  }

  const areaArcs: AreaArcs[] = areaPoints.map(d => {
    const arcs = pointsToArcs(d.points).map(arc => ({points: arc, line: lineEquation(arc)}))
    return {id: d.id, arcs: arcs}
  })

  const neighbors:[number,number][] = []
  
  possibleNeighborPolygons.forEach(pair => {
    const areaArcs1 = areaArcs[pair[0]]
    const areaArcs2 = areaArcs[pair[1]]
    if (areaArcsOverlap(areaArcs1, areaArcs2, 0)) {
      neighbors.push(pair)
    }
  })
  
  const neighborMap = new Map();
  [...neighbors.values()].forEach(d => {
    const pair = [areaArcs[d[0]].id,areaArcs[d[1]].id]
    if (pair[0] != pair[1]) {
      neighborMap.set(pair.join("_&_"), pair)
    }
    
  })

  return [...neighborMap.values()]
}

type LineSegment = [Position, Position]
function lineEquation(seg:LineSegment): LineEquation {
  const slope = (seg[1][1]-seg[0][1])/(seg[1][0]-seg[0][0])
  const intercept = seg[0][1] - slope*seg[0][0]
  return {m: slope, c: intercept}
}

function pointsBoundingBox(points: Position[]): [Position, Position] | null {
  const xExtent = d3.extent(points, d => d[0])
  const yExtent = d3.extent(points, d => d[1])
  if (xExtent[0] != null && xExtent[1] != null && yExtent[0] != null && yExtent[1] != null) {
    return [[xExtent[0], yExtent[0]],[xExtent[1],yExtent[1]]]
  } else {
    return null
  }
}

type Rectangle = [Position,Position]
function rectanglesOverlap(rect1:Rectangle, rect2:Rectangle, margin=0) {
  const rectWithMarg = [rect1[0].map(d => d-margin), rect1[1].map(d => d+margin)]
  if (rectWithMarg[1][0] < rect2[0][0] || rect2[1][0] < rectWithMarg[0][0]) return false 
  if (rectWithMarg[1][1] < rect2[0][1] || rect2[1][1] < rectWithMarg[0][1]) return false
  return true 
}

function pointsToArcs(points:Position[]) {
  const arcs: Arc[] = []
  for (let i = 0; i < points.length-1; i++) {
    arcs.push([points[i], points[i+1]])
  }
  return arcs 
}

function areaArcsOverlap(areaArcs1: AreaArcs, areaArcs2: AreaArcs, threshold=0) {
  for (const arc1 of areaArcs1.arcs) { 
    for (const arc2 of areaArcs2.arcs) {
      if (_segmentsOverlap(arc1.points, arc2.points, arc1.line, arc2.line) > threshold) {
        return true 
      }
    }
  }

  return false 
}

function _segmentsOverlap(seg1:LineSegment, seg2:LineSegment, 
    l1:LineEquation, l2:LineEquation, tolerance=0.00001) {
  const linesCollinear = Math.abs(l1.m - l2.m) < tolerance && Math.abs(l1.c - l2.c) < tolerance

  if (!linesCollinear) {
    return 0
  }

  const range = (d:number, min:number, max:number) => Math.max(min, Math.min(d, max))
  const relativeDistances = seg2.map(p => distanceAlongSegment(p, seg1))
  const trimmedDistances = relativeDistances.map(d => range(d, 0,1))
  const distance = Math.abs(trimmedDistances[0] - trimmedDistances[1])

  
  const length = Math.sqrt((seg1[0][0]-seg1[1][0])**2 + (seg1[0][1]-seg1[1][1])**2)
  return distance * length 
}

function distanceAlongSegment(p:Position, seg:LineSegment) {
  return (seg[0][0]-p[0])/(seg[0][0]-seg[1][0])
}