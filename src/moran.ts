import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm"
import {mean, deviation, sum} from "https://cdn.jsdelivr.net/npm/d3-array@3/+esm"

import type { FeatureCollection, } from "geojson"
import type { ID, Method } from "neighbors"
import {findNeighbors} from "./neighbors.js"

export {findNeighbors} from "./neighbors.js"


// TODO: Remove test code
// d3.json("/moranjs/data/dataWithMissing.geojson").then((result: any) => {
//   const featureCollection: FeatureCollection = result
//   const moranResult = calculateMoran(featureCollection, "adult_smoking").then(result => {
//     console.log(result)
//   })
//   // const printEvery = 50 / moranResult.localResults.length
//   // calculatePValues(moranResult, {permutations: 99, progressCallback: p => {
//   //   if (p % printEvery == 0) {
//   //     console.log(p)
//   //   }
//   //   return p
//   // }})
// })


// TODO: Remove debugging information
interface LocalMoranResult {
  id: string | number, value: number | undefined, z: number | undefined,
  lag?: number, moranLocal?: number,  p?: number, pCutoff?: number, label?: string,
  neighborZs?: (number | undefined)[]
}

interface MoranResult {
  localResults: LocalMoranResult[],
  moranGlobal: number, m2: number,
  weightMatrix: WeightMatrix,
  p?: number, pRef?: number[]
}


export async function calculateMoran(
  featureCollection: FeatureCollection, valueField: string, opts:{
    weightMatrix?: WeightMatrix, neighborMethod?: Method, localNormalize?: boolean,
    permutations?: undefined | number
  }={}) {

  Object.assign({
    neighborMethod: "Queen",
    localNormalize: true,
    permutations: undefined
  }, opts)

  if (opts.weightMatrix == null) {
    opts.weightMatrix = calculateWeightMatrix(featureCollection, opts.neighborMethod)
  }

  // Initialise basic results
  const localResults:LocalMoranResult[] = []
  featureCollection.features.forEach((feature, i) => {
    if (feature.properties != null) {
      const id = feature.id != null ? feature.id! : i

      let value = feature.properties[valueField]
      let z = undefined
      if (!opts.weightMatrix?.map.has(id)) {
        z = undefined
        value = undefined
      }

      localResults.push({
        id: id, value, z
      })
    }
  })

  // Calculate normalized values
  const valueMean = mean(localResults, d => d.value) as number
  const valueStd = deviation(localResults, d => d.value) as number
  
  localResults.forEach(localResult => {
    if (localResult.value != null) {
      localResult.z = (localResult.value - valueMean) / valueStd
    }
  })

  const localResultMap = new Map(localResults.map(d => [d.id, d]))

  // Calculate global and local moran values
  let moranGlobal = 0 
  let m2 = 0
  for (const localResult of localResults) {
    const areaWeightEntry = opts.weightMatrix?.getWeights(localResult.id)
    if (areaWeightEntry == undefined) continue

    const neighborIds = [...areaWeightEntry.keys()]
    const neighbors = [...neighborIds!].map(d => localResultMap.get(d))
    const neighborZs = neighbors.map(d => d!.z)
    const weights = [...areaWeightEntry!.values()]

    if (localResult.z != null) {
      const localMoranResult = unadjustedLocalMoran(localResult.z, neighborZs, weights)

      localResult.lag = localMoranResult.lag
      localResult.moranLocal = localMoranResult.moranLocal
      localResult.neighborZs = neighborZs

      if (localResult.moranLocal != undefined) {
        m2 += localResult.z**2
        moranGlobal += localResult.moranLocal 
      }
    }
    
  }
  moranGlobal = moranGlobal / m2

  if (opts.localNormalize) {
    for (const localResult of localResults) {
      if (localResult.moranLocal != undefined) {
        localResult.moranLocal = localResult.moranLocal / m2
      }
    }
  }

  const result = {
    localResults, moranGlobal, m2, weightMatrix: opts.weightMatrix
  }

  if (opts.permutations != undefined) {
    await calculatePValues(result, {permutations: opts.permutations})
  }

  return result 
}

export async function calculatePValues(moranResult:MoranResult, opts:{
  progressCallback?:(p:number)=>number,
  permutations?: number,
}={}) {
  Object.assign(opts, {
    progressCallback: (p:number) => p,
    permutations: 999,
  })

  const permutations = opts.permutations !
  const progressCallback = opts.progressCallback!

  // This is approximating proper permutation
  const localResults = moranResult.localResults
  const permutes:LocalMoranResult[][] = []
  for (let i = 0; i < permutations; i++) {
    const permute = [...localResults]
    d3.shuffle(permute)
    permutes.push(permute)
  }

  const maxNeighbors = d3.max([...moranResult.weightMatrix.map.values()], d => d.size)
  if (maxNeighbors! * 2 > moranResult.localResults.length) {
    // TODO: Fix this, don't need to throw error.
    new Error("Max neighbors * 2 greater than number of locations")
  }

  
  localResults.forEach((localResult,i) => {
    if (localResult.z != undefined) {
      // const zValues = localResults.map(d => d.z)
      // zValues.splice(i, 1)
      const weightRow = moranResult.weightMatrix.getWeights(localResult.id)
      
      let moreExtremeCount = 0

      const bsLocalMorans = []
      if (weightRow != undefined) {
        const weights = [...weightRow.values()]
        for (let j = 0; j < permutations; j++) {
          const permute = permutes[j]
          let sample = permute.slice(0, weights.length)
          if (sample.some(d => d.id == localResult.id)) {
            sample = permute.slice(weights.length, weights.length+weights.length)
          }
          const neighborZs = sample.map(d => d.z)
          const moranLocal = unadjustedLocalMoran(localResult.z, neighborZs, weights)
          if (moranLocal.moranLocal != undefined && localResult.moranLocal != undefined) {
            bsLocalMorans.push(moranLocal.moranLocal)
            if (Math.abs(moranLocal.moranLocal) >= Math.abs(localResult.moranLocal)) {
              moreExtremeCount++
            }
          }

          
        }
      }

  
      const actualLocalMoran = localResult.moranLocal
      if (actualLocalMoran != undefined) {
        const unadjustedLocalMoran = actualLocalMoran// * moranResult.m2
        // const refLocalMorans = bsLocalMorans.filter(
        //   actualLocalMoran >= 0 ? d => d > 0 : d => d < 0).map(d => Math.abs(d))
        // refLocalMorans.sort((a,b) => b - a)
        // const minIndex = refLocalMorans.findIndex(d => Math.abs(actualLocalMoran) > d)
        //localResult.p = (minIndex + 1) / (permutations + 1)
        
        let p = pByReference(unadjustedLocalMoran, bsLocalMorans, permutations)

      //   // TODO: Remove
      //  // if (localResult.id == "19077") {
      //     //console.log(unadjustedLocalMoran, bsLocalMorans, permutations, p)

      //     const ref = bsLocalMorans
      //     const val = unadjustedLocalMoran
      //     //const permutations = per

      //     let referenceArr = ref.filter(d => Math.sign(d) == Math.sign(val))
      //     referenceArr = referenceArr.map(d => Math.abs(d))
      //     referenceArr.sort((a,b) => a - b)
      //     const absVal = Math.abs(val)
      //     const minIndex = referenceArr.findIndex(d => absVal < d)

      //     //console.log(referenceArr, absVal, minIndex, (ref.length- minIndex + 1) /  (permutations + 1) )

      //     p = (ref.length - minIndex + 1) /  (ref.length + 1)
      //   //}

      //   // TODO: Consider...
      //   p = (moreExtremeCount+1)/(permutations+1)/2

        // const ref = bsLocalMorans
        // const val = unadjustedLocalMoran

        // let referenceArr = ref.filter(d => Math.sign(d) == Math.sign(val))
        // referenceArr = referenceArr.map(d => Math.abs(d))
        // referenceArr.sort((a,b) => a - b)
        // const absVal = Math.abs(val)
        // let minIndex = referenceArr.findIndex(d => absVal < d)
        // if (minIndex < 0) {
        //   minIndex = ref.length
        // }

        // p = (referenceArr.length - minIndex + 1) /  (permutations + 1)


        // if (localResult.id == "17043") {
        //  // console.log(absVal, referenceArr, minIndex, p )
        //  console.log(minIndex, (referenceArr.length - minIndex + 1) ,  (permutations + 1), p )
        // }


        localResult.p = p
        if (p != null) {
          localResult.pCutoff = [0.0001, 0.001, 0.01, 0.05].find(d => p < d)
          if (p < 0.05) {
              let label = ""
              label = label + (localResult.z >= 0 ? "High" : "Low")
              label = label + (localResult.lag! >= 0 ? "-High" : "-Low")
              localResult.label = label
            } else {
              localResult.label = "Not significant"
            }
        }
        
      }
    }
    progressCallback(i/localResults.length)
  })

  const bsGlobalMorans = []
  const ids = localResults.map(d => d.id)
  const zs = localResults.map(d => d.z)
  for (let i = 0; i < permutations; i++) {
    d3.shuffle(zs)
    const zMap = new Map(ids.map((d,i) => [d,zs[i]]))
    let globalMoran = 0
    for (const id of ids) {
      const weights = moranResult.weightMatrix.getWeights(id)
      const z = zMap.get(id)
      if (weights != null && z != null) {
        const neighborIds = [...weights.keys()]
        const weightValues = [...weights.values()]
        const neighborZs = neighborIds.map(id => zMap.get(id))
        const moranLocalResult = unadjustedLocalMoran(z, neighborZs, weightValues)
        if (moranLocalResult.moranLocal != undefined) {
          globalMoran += moranLocalResult.moranLocal
        }
      }
    }
    bsGlobalMorans.push(globalMoran)
  }

  const unadjustedGlobalMoran = moranResult.moranGlobal*moranResult.m2
  const globalP = pByReference(unadjustedGlobalMoran, bsGlobalMorans, permutations)
  moranResult.pRef = bsGlobalMorans

  moranResult.p = globalP
  
  progressCallback(1)
  return moranResult
}

// I am not certain if this function is correct, but it gives results close to GeoDa's output. The calculation should 
// be straightforward, but GeoDa's results are exactly half what I would expect (and no p-values are >0.5). I reckon
// they are folding the p-values to [0,0.5], but I don't quite understand why this is justified given that Moran's I
// is essentially a 2-sided test. I think in only comparing the values to ones of the same sign, we are essentially 
// doing a 1-sided test by design, whereas GeoDa is forcing it? I am not sure.
function pByReference(val:number, ref:number[], permutations:number) {
  // TODO: Replace this with just a count of more extreme, if that can give the same results. 

  // let referenceArr = ref.filter(d => Math.sign(d) == Math.sign(val))
  // referenceArr = referenceArr.map(d => Math.abs(d))
  // referenceArr.sort((a,b) => b - a)
  // const minIndex = referenceArr.findIndex(d => Math.abs(val) > d)
  // return (minIndex + 1) /  (permutations + 1)

  let referenceArr = ref.filter(d => Math.sign(d) == Math.sign(val))
  referenceArr = referenceArr.map(d => Math.abs(d))
  referenceArr.sort((a,b) => a - b)
  const absVal = Math.abs(val)
  let minIndex = referenceArr.findIndex(d => absVal < d)
  if (minIndex < 0) {
    minIndex = ref.length
  }
  return (referenceArr.length - minIndex + 1) /  (permutations + 1)
}

// TODO: Undefined behavior if no valid neighbors
function unadjustedLocalMoran(z: number, neighborZs: (number|undefined)[], weights: number[]) {
  const validNeighborZs:number[] = []
  let validWeights:number[] = []

  for (let i = 0; i < neighborZs.length; i++) {
    const neighborZ = neighborZs[i]
    if (neighborZ != null && Number.isFinite(neighborZ)) {
      validNeighborZs.push(neighborZ)
      validWeights.push(weights[i])
    }
  }

  const weightSum = sum(validWeights)
  validWeights = validWeights.map(w => w / weightSum)

  if (validWeights.length > 0) {
    let lag = 0
    for (let i = 0; i < validWeights.length; i++) {
      lag += validNeighborZs[i] * validWeights[i]
    }
    return {moranLocal: z * lag, lag: lag}
  } else {
    return {moranLocal: undefined, lag: undefined}
  }
  
}

export function calculateWeightMatrix(featureCollection:FeatureCollection, method:Method = "queen") {
  const neighbors = findNeighbors(featureCollection, method)
  const neighborPairMap = d3.group(neighbors, d => d[0])
  const weightMatrix = new WeightMatrix()

  for (const pairs of neighborPairMap.values()) {
    pairs.forEach(pair => weightMatrix.set(pair[0], pair[1], 1/pairs.length))
  }

  return weightMatrix
}

export class WeightMatrix {
  map: Map<ID, Map<ID, number>>

  constructor(map?:Map<ID, Map<ID, number>>) {
    if (map) {
      this.map = map
    } else {
      this.map = new Map()
    }
  }

  getWeights(id: ID) {
    return this.map.get(id)
  }

  getWeight(id1:ID, id2: ID) {
    return this.map.get(id1)?.get(id2)
  }

  set(id1:ID, id2: ID, value: number) {
    let mapRow = this.map.get(id1)
    if (mapRow == undefined) {
      mapRow = new Map()
      this.map.set(id1, mapRow)
    }
    mapRow.set(id2, value)
  }
}