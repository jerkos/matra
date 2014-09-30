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

class MassTrace {
	double mz;
	Peak apexPeak;
	DList!Peak peaks;
	size_t size;

	this() {
		size=0;
		apexPeak=null;
	}

	nothrow void append(Peak p) {
		peaks.insertFront(p);
		++size;
	}

	nothrow void appendLeft(Peak p) {
		peaks.insertBack(p);
		++size;
	}

	nothrow size_t length() {
		return size;
	}
}

class PeakIndex {
	static double BIN_SIZE = 0.25;
	static double INV_BIN_SIZE = 4;

	//Peak[] peaks;
	Peak[][int] index; 
	double minMz;
	double maxMz;

	this() { //Peak[] peaks_) {
		//peaks = peaks_;
		minMz = 60;
		maxMz = 895;

		//update(peaks_);
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
    	writeln("length rtByScanId:", rtByScanId.length);
    	sqlite3_finalize(stmt);
    	sqlite3_close(db);
	}

	void extract(out MassTrace[] massTraces) {
		auto rsIt = new test.RunSliceIterator(filename, 1);

		Peak[] allPeaks;
		//allPeaks.reserve(10000000);
		PeakIndex[int] peaksIdxByScanId;
		
		foreach(int scanId; scanIds)
			peaksIdxByScanId[scanId] = new PeakIndex();
		

		auto currTuple = rsIt.next();
		test.RSHeader currRsh = currTuple[0];
		test.ScanSlice[] currSs = currTuple[1];

		test.RSHeader prevRsh = null;
		test.ScanSlice[] prevSs;

		test.RSHeader nextRsh = null;
		test.ScanSlice[] nextSs;

		foreach(ref test.ScanSlice s; currSs) {
			auto peaks = s.toPeaks(rtByScanId);
			peaksIdxByScanId[s.scanId].update(peaks);
			allPeaks ~= peaks;
		}

		int[Peak] alreadyUsedPeaks;

		while (rsIt.hasNext()) {
			auto nextTuple = rsIt.next();
			nextSs = nextTuple[1];
			nextRsh = nextTuple[0];

			writeln("Processing runslice #", nextRsh.id);

			foreach(ref test.ScanSlice s; nextSs)
				allPeaks ~= s.toPeaks(rtByScanId);

			Peak[] notSeenYet = allPeaks.filter!(x=> x !in alreadyUsedPeaks)
										.array
										.sort!((a, b)=> a.intensity > b.intensity)
										.array;


			foreach(ref Peak peak; notSeenYet) {
				if (peak in alreadyUsedPeaks)
					continue;
			
				MassTrace massTrace = new MassTrace();
				massTrace.apexPeak = peak;
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
				if (massTrace.length() > 2)
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

			prevRsh = currRsh;
			prevSs = currSs;

			currRsh = nextRsh;
			currSs = nextSs;

		}
	}
	//}
}

void main() {
	auto file = "C:\\Users\\Marco\\Desktop\\X20140626_006DP_pos_122.raw.mzdb";
	auto rsIt = new test.RunSliceIterator(file, 1);

	auto currentTime = Clock.currTime();
	MassTrace[] massTraces;
	MassTraceExtractor mtd = new MassTraceExtractor(file, 15.0, 1);
	mtd.extract(massTraces);
	writeln("elapsed time:", Clock.currTime() - currentTime);
}