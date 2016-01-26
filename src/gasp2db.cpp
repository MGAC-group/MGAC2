#include "gasp2db.hpp"


using namespace std;


GASP2db::GASP2db(string name) {

	dbconn=NULL;
	path = name;
	openState=false;

}

int GASP2db::connect() {
	int state;
	if(!openState) {
		state = sqlite3_open(path.c_str(), &dbconn);
		if(state) {
			//something bad happened on opening the SQL
			//time to panic!
			cout << mark() << "ERROR: SQLITE RETURN WITH ERROR CODE " << state << ". ABORTING!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
			exit(1);
		}
		else {
			openState=true;
		}
	}

	return 1;
}

int GASP2db::disconnect() {

	if(openState) {
		//if close returns as SQLITE_BUSY, it's okay
		//we just try to close later
		if(sqlite3_close(dbconn) == SQLITE_OK) {
			openState=false;
			return 1;
		}
		else {
			//openState=true;
			return 0;
		}
	}
	else
		return 1;

}

int GASP2db::load(string name) {
	if(!openState) {
		path=name;

		return 1;
	}
	else
		return 0;
}

void GASP2db::init() {

	char * err = 0;

	if(connect()) {

		sqlite3_exec(dbconn, "BEGIN TRANSACTION", NULL, NULL, &err);

		//the I prefix indicates initial, otherwise the value is current
		sqlite3_exec(dbconn,
				"CREATE TABLE IF NOT EXISTS structs( "
				"id TEXT PRIMARY KEY,"
				"parentA TEXT, parentB TEXT, generation INT, version INT,"
				"energy REAL, force REAL, pressure REAL, spacegroup INT, cluster INT,"
				"axis INT, type INT, centering REAL, subtype REAL,"
				"state INT, error INT, time INT, steps INT,"
				"a REAL, b REAL, c REAL, al REAL, bt REAL, gm REAL, ra REAL, rb REAL, rc REAL,"
				"Ia REAL, Ib REAL, Ic REAL, Ial REAL, Ibt REAL, Igm REAL, Ira REAL, Irb REAL, Irc REAL,"
				"tA REAL, tB REAL, tC REAL, mB REAL, rhmC REAL,"
				"ItA REAL, ItB REAL, ItC REAL, ImB REAL, IrhmC REAL,"
				"coord BLOB, Icoord BLOB"
				")", NULL, NULL, &err);

		sqlite3_exec(dbconn, "CREATE TABLE IF NOT EXISTS types(id INT PRIMARY KEY, name TEXT, xml TEXT)", NULL, NULL, &err);

		sqlite3_exec(dbconn, "CREATE TABLE IF NOT EXISTS inputs(id INT PRIMARY KEY, xml TEXT)", NULL, NULL, &err);

		sqlite3_exec(dbconn, "END TRANSACTION", NULL, NULL, &err);

		disconnect();
	}

}


//sets the INITIAL RECORDS of the structs table
bool GASP2db::create(GASP2pop pop) {

	char * err = 0;
	sqlite3_stmt * stm;

	if(connect()) {

		sqlite3_prepare_v2(dbconn, "INSERT INTO structs ("
				"id, parentA, parentB, generation, version,"
				"energy, force, pressure, spacegroup, cluster,"
				"axis, type, centering, subtype,"
				"state, error, time, steps,"
				"Ia, Ib, Ic, Ial, Ibt, Igm, Ira, Irb, Irc,"
				"ItA, ItB, ItC, ImB, IrhmC,"
				"Icoord"
				") VALUES ("
				"@id, @pa, @pb, @gen, @ver,"
				"@en, @fr, @pr, @spcg, @cl,"
				"@ax, @typ, @cen, @sub,"
				"@st, @er, @tm, @steps,"
				"@Ia, @Ib, @Ic, @Ial, @Ibt, @Igm, @Ira, @Irb, @Irc,"
				"@ItA, @ItB, @ItC, @ImB, @IrhmC, "
				"@Icoord"
				")", 1, &stm, NULL);

		sqlite3_exec(dbconn, "BEGIN TRANSACTION", NULL, NULL, &err);

		for(int i = 0; i < pop.size(); i++)
			pop.indv(i)->sqlbindCreate(stm);

		sqlite3_exec(dbconn, "END TRANSACTION", NULL, NULL, &err);

		disconnect();

		return true;

	}

	return false;


}



















