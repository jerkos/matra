module extractor;

import std.algorithm;
import std.datetime;
import std.stdio;
import std.string;
import std.container;
import std.math;
import std.array;
import etc.c.sqlite3;
import std.conv;

import model;

/++
 + Mass trace class 
 +/
class MassTrace {
	Peak apexPeak;
	DList!Peak peaks;
	size_t size;

	this() {
		size=0;
		apexPeak=null;
	}


	nothrow void setApexPeak(ref Peak p) {
		apexPeak = p;
		++size;
	}

	nothrow void append(ref Peak p) {
		peaks.insertFront(p);
		++size;
	}

	nothrow void appendLeft(ref Peak p) {
		peaks.insertBack(p);
		++size;
	}

}

/++
 + 
 +/
class PeakIndex {
	static double BIN_SIZE = 0.25;
	static double INV_BIN_SIZE = 4;

	Peak[][int] index; 
	double minMz;
	double maxMz;

	this() { 
		minMz = 60.0;
		maxMz = 895.0;
	}

	int empty() {
		return index.length;
	}

	void update(Peak[] peaks_) {
		foreach(ref Peak p; peaks_) {
			auto bin = to!int(p.mz * PeakIndex.INV_BIN_SIZE);
			index[bin] ~= p;
		}
	}

	void removeByMass(in double minMz, in double maxMz) {
		auto minIdx = to!int(minMz * INV_BIN_SIZE);
		auto maxIdx = to!int(maxMz * INV_BIN_SIZE);

		for(int i = minIdx; i <= maxIdx; ++i) {
			if (i in index)
				index.remove(i);
		}
	}

	Peak getNearestPeak(in double mzVal, in double tol, ref int[Peak] alreadyUsed) {
		
		auto bin = to!int(mzVal * INV_BIN_SIZE);

		Peak[] matchingPeaks;
		for (int i = bin - 1; i <= bin + 1; ++i) {
			if (i in index) {
				matchingPeaks ~= index[i];
			}
		}
		
		if (matchingPeaks.empty())
			return null;

		Peak nearestPeak = null;
		double minDist = 1e6;
		foreach(ref Peak p; matchingPeaks) {
			if (p in alreadyUsed)
				continue;
			auto diff = abs(p.mz - mzVal);
			if (diff > tol)
				continue;
			if (diff < minDist) {
				nearestPeak = p;
				minDist = diff;
			}
		}
		return nearestPeak;
	}
}

/++
 +
 +/
class MassTraceExtractor {
	double mzTolPPM;
	int gapAllowed;
	float[int] rtByScanId;
	string filename;
	int[] scanIds;

	this(string _filename, double _mzTolPPM, int _gapAllowed) {
		filename = _filename;
		mzTolPPM = _mzTolPPM;
		gapAllowed = _gapAllowed;
		fillRtByScanId();
		scanIds = rtByScanId.keys.sort;
	}

	void fillRtByScanId() {
		sqlite3* db;
		sqlite3_stmt* stmt;

		sqlite3_open_v2(toStringz(filename), &db, SQLITE_OPEN_READONLY, null);
		
		auto sqlQuery = "SELECT id, time FROM spectrum WHERE ms_level = 1;";
		sqlite3_prepare_v2(db, toStringz(sqlQuery), -1, &stmt, null);
    	while(sqlite3_step(stmt) == SQLITE_ROW) {
    		int id = sqlite3_column_int(stmt, 0);
    		float rt = sqlite3_column_double(stmt, 1);
    		rtByScanId[id] = rt;
    	}
    	writefln("length rtByScanId: %d", rtByScanId.length);
    	sqlite3_finalize(stmt);
    	sqlite3_close(db);
	}


	void update(ref ScanSlice[] currSs, 
				ref PeakIndex[int] peaksIdxByScanId, 
				ref Peak[] allPeaks) {
		
		foreach(ref ScanSlice s; currSs) {
			auto peaks = s.toPeaks(rtByScanId);
			peaksIdxByScanId[s.scanId].update(peaks);
			allPeaks ~= peaks;
		}
	}

	/++
	 + main function called here
	 +/
	void extract(out MassTrace[] massTraces) {
		auto rsIt = new RunSliceIterator(filename, 1);

		/// maybe could reserve this one ?
		///allPeaks.reserve(10000000);
		Peak[] allPeaks;
		
		/// setting peak index for each scanID
		PeakIndex[int] peaksIdxByScanId;
		foreach(ref int scanId; scanIds)
			peaksIdxByScanId[scanId] = new PeakIndex();
		
		/// first iteration setting to currRunSlice
		auto currTuple = rsIt.next();
		RSHeader currRsh = currTuple[0];
		ScanSlice[] currSs = currTuple[1];

		/// updating peaksByScanID and allpeaks
		//foreach(ref ScanSlice s; currSs) {
		//	auto peaks = s.toPeaks(rtByScanId);
		//	peaksIdxByScanId[s.scanId].update(peaks);
		//	allPeaks ~= peaks;
		//}
		update(currSs, peaksIdxByScanId, allPeaks);

		/// setting next runslice as we are iterating on a slide window
		RSHeader prevRsh = null;
		ScanSlice[] prevSs;

		RSHeader nextRsh = null;
		ScanSlice[] nextSs;

		/// hashmap mimic set python datastucture
		int[Peak] alreadyUsedPeaks;

		/// starting main while loop
		while (rsIt.hasNext()) {

			/// fetch next
			auto nextTuple = rsIt.next();
			nextSs = nextTuple[1];
			nextRsh = nextTuple[0];

			writeln("Processing runslice #", nextRsh.id);

			update(nextSs, peaksIdxByScanId, allPeaks);


			Peak[] notSeenYet = allPeaks.filter!(x=> x !in alreadyUsedPeaks)
										.array
										.sort!((a, b)=> a.intensity > b.intensity)
										.array;


			foreach(ref Peak peak; notSeenYet) {
				if (peak in alreadyUsedPeaks)
					continue;
			
				MassTrace massTrace = new MassTrace();
				massTrace.setApexPeak(peak);
				alreadyUsedPeaks[peak] = 0;
				
				auto scanId = peak.scanId;

				auto currMz = peak.mz;

				auto tol = currMz * mzTolPPM * 15;

				auto scanIdIdx = countUntil(scanIds, scanId);

				auto prevScanIdIdx = scanIdIdx - 1;
				int prevScanId = -1;
				if (prevScanIdIdx >= 0) {
					prevScanId = scanIds[prevScanIdIdx];
				}

				auto nextScanIdIdx = scanIdIdx + 1;
				
				int nextScanId = -1;
				if (nextScanIdIdx < scanIds.length) {
					nextScanId = scanIds[nextScanIdIdx];
				}

				int rightGap = 0;
				while (1) { //nextScanId != -1) {
					if (nextScanId !in peaksIdxByScanId)
						break;
					auto peakIdx = peaksIdxByScanId[nextScanId];

					auto p = peakIdx.getNearestPeak(currMz, tol, alreadyUsedPeaks);
					//Peak p;
					if (p is null) {
						rightGap++;
						if (rightGap > gapAllowed)
							break;
					} else {
						massTrace.append(p);
						alreadyUsedPeaks[p] = 0;
					}
					nextScanIdIdx++;
					if (nextScanIdIdx < scanIds.length) {
						nextScanId = scanIds[nextScanIdIdx];
					} else {
						break;
					}
				}

				int leftGap = 0;
				while (1) {
					if (prevScanId !in peaksIdxByScanId)
						break;
					auto peakIdx = &peaksIdxByScanId[prevScanId];
					auto p = peakIdx.getNearestPeak(currMz, tol, alreadyUsedPeaks);
					//Peak p;
					if (p is null) {
						leftGap ++;
						if (leftGap > gapAllowed)
							break;
					} else {
						massTrace.appendLeft(p);
						alreadyUsedPeaks[p] = 0;
					}

					prevScanIdIdx--;
					if (prevScanIdIdx >= 0) {
						prevScanId = scanIds[prevScanIdIdx];
					} else {
						break;
					}
				}
		
				if (massTrace.size > 5)
					massTraces ~= massTrace;
			}

			allPeaks = allPeaks.filter!(x=> prevRsh !is null ? 
				x.rsId == prevRsh.id : false).array; 

			// remove peaks from index
			if (prevRsh !is null) {
				auto maxmz = 60 + (prevRsh.id * 5) - 1;
				auto minmz = maxmz - 5;
				foreach(ref PeakIndex pi; peaksIdxByScanId.byValue()) {
					pi.removeByMass(minmz, maxmz);
				}
			}

			/// Updating stuffs
			prevRsh = currRsh;
			prevSs = currSs;

			currRsh = nextRsh;
			currSs = nextSs;

		}
	}
	//}
}

void main() {
	auto file = "D:\\Utilisateurs\\Marc\\Desktop\\developpement\\082DBB.raw.mzDB";
	auto rsIt = new RunSliceIterator(file, 1);

	auto currentTime = Clock.currTime();
	MassTrace[] massTraces;
	MassTraceExtractor mtd = new MassTraceExtractor(file, 15.0, 1);
	mtd.extract(massTraces);
	writeln("mass traces extracted:", massTraces.length);
	writeln("elapsed time:", Clock.currTime() - currentTime);
}