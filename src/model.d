module model;
pragma(lib, "lib/sqlite3");

import std.stdio;
import std.conv;
import std.string;
import std.datetime;
import std.typecons;
import etc.c.sqlite3;

class Peak {
	double mz;
	float intensity;
	float rt;
	int rsId;
	int scanId;

	this(double mz_, float intensity_, float rt_, int rsId_, int scanId_) {
		mz = mz_;
		intensity = intensity_;
		rt = rt_;
		rsId = rsId_;
		scanId = scanId_;
	}
}

class BoundingBox {
	int id;
	ubyte* data;
	int blobLength;
	int runSliceId;
	int firstSpectrumId;
	int lastSpectrumId;

	this(int id_, ubyte* data_, int blobLength_, int runSliceId_, int firstSpectrumId_, int lastSpectrumId_) {
		id = id_;
		data = data_;
		blobLength = blobLength_;
		runSliceId = runSliceId_;
		firstSpectrumId = firstSpectrumId_;
		lastSpectrumId = lastSpectrumId_;
	}
}

class RSHeader {
	int id;
	int msLevel;
	int number;
	float beginMz;
	float endMz;

	this(int id_, int msLevel_, int number_, float beginMz_, float endMz_) {
		id = id_;
		msLevel = msLevel_;
		number = number_;
		beginMz = beginMz_;
		endMz = endMz_;
	}
}

class ScanSlice {
	int scanId;
	int rsId;
	double[] mzs;
	float[] ints;

	this(int scanId_, int rsId_, ref double[] mzs_, ref float[] ints_) {
		scanId = scanId_;
		rsId = rsId_;
		mzs = mzs_;
		ints = ints_;
	}

	Peak[] toPeaks(ref float[int] rtByScanId) {
		Peak[] peaks;
		peaks.length = mzs.length;
		for (int i = 0; i < mzs.length; ++i) {
			peaks[i] = new Peak(mzs[i], ints[i], rtByScanId[scanId], rsId, scanId);
		}
		return peaks;
	}
}

class RunSliceIterator {
	static string sqlQuery = "SELECT bounding_box.* FROM bounding_box,
                   run_slice WHERE run_slice.ms_level = 1
                   AND bounding_box.run_slice_id = run_slice.id
                   ORDER BY run_slice.begin_mz;";

    static string headerSqlQuery = "select id, ms_level, number, 
    begin_mz, end_mz from run_slice where ms_level = 1";

    sqlite3* db;
    sqlite3_stmt* stmt;
    int msLevel;
    BoundingBox firstRSBB;
    ScanSlice[] currScanSlices;
    RSHeader[int] rsHeaderById;

    this(string filename, int msLevel) {
    	msLevel = msLevel;
    	sqlite3_open_v2(toStringz(filename), &db, SQLITE_OPEN_READONLY, null);
    	
    	sqlite3_exec(db, toStringz("PRAGMA synchronous=OFF;
                         PRAGMA journal_mode=OFF;
                         PRAGMA temp_store=3;
                         PRAGMA cache_size=8000;
                         PRAGMA page_size=4096;"), null, null, null);

    	fillRSHeaderById();
    	
    	sqlite3_prepare_v2(db, toStringz(sqlQuery), -1, &stmt, null);
    	sqlite3_step(stmt);

    	int bbId = sqlite3_column_int(stmt, 0);
    	ubyte* data = cast(ubyte*)sqlite3_column_blob(stmt, 1);
    	int rsId = sqlite3_column_int(stmt, 2);
    	int firstSpecId = sqlite3_column_int(stmt, 3);
    	int lastSpecId = sqlite3_column_int(stmt, 4);

    	int blobLength = sqlite3_column_bytes(stmt, 1);
    	firstRSBB = new BoundingBox(bbId, data, blobLength, rsId, firstSpecId, lastSpecId);
    	
    	//sqlite3_finalize(stmt);
    }

    ~this() {
    	sqlite3_close(db);
    }

    void fillRSHeaderById() {
    	sqlite3_prepare_v2(db, toStringz(headerSqlQuery), -1, &stmt, null);
		while (sqlite3_step(stmt) != SQLITE_DONE) {
			int id = sqlite3_column_int(stmt, 0);
			int msLevel_ = sqlite3_column_int(stmt, 1);
			int number = sqlite3_column_int(stmt, 2);
			float beginMz = sqlite3_column_double(stmt, 3);
			float endMz = sqlite3_column_double(stmt, 4);
			rsHeaderById[id] = new RSHeader(id, msLevel_, number, beginMz, endMz); 
		}
		sqlite3_finalize(stmt);    
    }


    void initIter() {
    	currScanSlices = bb2ScanSlice(firstRSBB);
    	while (1) {
    		int r = sqlite3_step(stmt);
    		if (r == SQLITE_DONE) {
    			firstRSBB = null;
    			break;
    		}

    		int bbId = sqlite3_column_int(stmt, 0);
	    	ubyte* data = cast(ubyte*)sqlite3_column_blob(stmt, 1);
	    	int rsId = sqlite3_column_int(stmt, 2);
	    	int firstSpecId = sqlite3_column_int(stmt, 3);
	    	int lastSpecId = sqlite3_column_int(stmt, 4);

	    	int blobLength = sqlite3_column_bytes(stmt, 1);
	    	auto bb = new BoundingBox(bbId, data, blobLength, rsId, firstSpecId, lastSpecId);

	    	//sqlite3_finalize(stmt);

	    	if (bb.runSliceId == firstRSBB.runSliceId) {
	    		currScanSlices ~= bb2ScanSlice(bb);
	    	} else {
	    		firstRSBB = bb;
	    		break;
	    	}
    	}
    }

    bool hasNext() {    	
    	if (firstRSBB is null) {
    		sqlite3_finalize(stmt);
    		return false;
    	}
    	return true;
    }

    auto next() {
    	initIter();
    	auto rsId = currScanSlices[0].rsId;
    	return tuple(rsHeaderById[rsId], currScanSlices);
    }

    ScanSlice[] bb2ScanSlice(ref BoundingBox bb) {
    	int offset = 0;
    	ScanSlice[] ss;
    	while (offset < bb.blobLength) {
    		int scanId = *cast(int*)bb.data[offset .. offset + 4];
    		offset += 4;
    		int nbPeaks = *cast(int*)bb.data[offset .. offset + 4];
    		offset += 4;
    		if (!nbPeaks) {
    			continue;
    		}

    		double[] mzs; mzs.length = nbPeaks;
    		float[] ints_; ints_.length = nbPeaks;
    		for (int i=0; i < nbPeaks; ++i) {
    			double mz = *cast(double*)bb.data[offset .. offset + 8];
    			offset += 8;
    			float ints = *cast(float*)bb.data[offset .. offset + 4];
    			offset += 4;
    			//skip widths
    			offset += 8;
    			mzs[i] = mz;
    			ints_[i] = ints;
    		}	
    		ss ~= new ScanSlice(scanId, bb.runSliceId, mzs, ints_);
    	}
    	return ss;
    }

}

/*void main() {
	auto file = "C:\\Users\\Marco\\Desktop\\X20140626_006DP_pos_122.raw.mzdb";
	auto rsIt = new RunSliceIterator(file, 1);
	//writeln("size:", rsIt.rsHeaderById.length);
	auto currentTime = Clock.currTime();
	int c = 0;
	while (rsIt.hasNext()) {
		ScanSlice[] ss = rsIt.next()[1];
		foreach (ref ScanSlice s; ss) {
			c += s.mzs.length;
		}
	}
	writeln("data points count:", c);
	writeln("elapsed time:", Clock.currTime() - currentTime);
}*/