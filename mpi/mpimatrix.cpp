/**\file
Реализация пересылки подматриц.
При включённом useBigSends происходит сериализация матрицы, и пересылка её в сериализованном виде.
При выключенном - пересылка по отдельным строкам.
Понятия попарной и широковещательных персылок представляются в виде классов PeerConnector и BCastConnector,
предоставляющих одинаковый интерфейс в виде методов send и recv.
*/


#include "types.h"
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <cstdlib>
#include <vector>
using namespace std;
namespace F4MPI{

///число переменных типа int, необходимое для помещения данных RowElement
static const unsigned INTS_IN_ROW_ELEMENT=1+(sizeof(RowElement)-1)/sizeof(int);
///тип, представляющий сериализованные данные (матрицу)
typedef vector<int> SerialData;

/**сериализация матрицы
Записывает набор строк [\a from, \a to) в сериализованном виде в данные \a data.
Формат представления состоит из последовательно размещённых перемееных типа int:
\arg число строк в матрице
\arg число ненулевых элементов для каждой строки
\arg последовательно размещённые данные строки
Исключение представляет пустая матрица (0 строк), которая представляется в виде данных нулевой длинны
*/
void serializeSubMatrix(SerialData& data, CMatrix::const_iterator from, CMatrix::const_iterator to){
	if (to==from){//Пустая матрица имеет особое представление - пустой вектор
		data.clear();
		return;
	}
	size_t totalSize=1;//место под число строк
	CMatrix::const_iterator cur;
	for (cur=from;cur!=to;++cur){
		totalSize+=1+INTS_IN_ROW_ELEMENT*(cur->size());//место под число элементов строки и сами элементы
	}
	data.resize(totalSize);
	int lastFree=0;
	data[lastFree]=to-from;
	++lastFree;
	for (cur=from;cur!=to;++cur,++lastFree){
		data[lastFree]=cur->size();//число элементов строк
	}
	for (cur=from;cur!=to;++cur){
		memcpy(&data[lastFree],&cur->front(),cur->size()*sizeof(cur->front()));//данные строки
		lastFree+=INTS_IN_ROW_ELEMENT*cur->size();
	}
}

/**десериализация с сохранением по заданному итератору
Записывает набор строк представленный данными \a data в строки, начиная со строки, на которую указывает итератор \a i.
*/
void deSerializeSubMatrix(const SerialData& data, CMatrix::iterator i){
	if (data.empty()) return;
	int readPos=0;
	int numRows=data[readPos];
	++readPos;
	int lastReadPos=readPos+numRows;
	CMatrix::iterator j;
	for (j=i;readPos!=lastReadPos;++readPos,++j){
		j->resize(data[readPos]);
	}
	lastReadPos=data.size();
	for (j=i;readPos!=lastReadPos;++j){
		memcpy(&j->front(),&data[readPos],j->size()*sizeof(j->front()));//данные строки
		readPos+=INTS_IN_ROW_ELEMENT*j->size();
	}
}

/**десериализация в матрицу
Добавляет набор строк представленный данными \a data в матрицу \a m.
*/
void deSerializeToMatrix(const SerialData& data, CMatrix& m){
	if (data.empty()){
		return;
	}
	m.resize(m.size()+data[0]);
	deSerializeSubMatrix(data,m.end()-data[0]);
}

///Реализует попарные коммуникации MPI
struct PeerConnector{
	PeerConnector(int peer):peerRank(peer){}
	int peerRank;///<ранг процесса в MPI_COMM_WORLD с которым идёт обмен 
	///отправляет данные на которые указывает \a ptr, считая что там \a size элементов типа \a dataType
	void send(void* ptr, size_t size, MPI_Datatype dataType){
		//MPI_процедуры принимают в качечтсве размера данных int(?).
		//Возможно, здесь следует проверять это ограничение м по необходимости резать данные
		MPI_Send(ptr, size, dataType, peerRank, 0, MPI_COMM_WORLD);
	}
	///получает и записывает по указателю \a ptr, данные из \a size элементов типа \a dataType
	void recv(void* ptr, size_t size, MPI_Datatype dataType){
		MPI_Status dummy;
		MPI_Recv(ptr, size, dataType, peerRank, 0, MPI_COMM_WORLD, &dummy);
	}
};

///Реализует коллективные коммуникации MPI
struct BCastConnector{
	BCastConnector(int sender, MPI_Comm commToUse):senderRank(sender),comm(commToUse){}
	int senderRank;///<ранг процесса, отправляющего данные при коллективном обмене
	MPI_Comm comm;///<Коммуникатор, определяющий набор процессов для коллективного обмена
	///отправляет данные на которые указывает \a ptr, считая что там \a size элементов типа \a dataType
	void send(void* ptr, size_t size, MPI_Datatype dataType){
		MPI_Bcast(ptr, size, dataType, senderRank, comm);
	}
	///получает и записывает по указателю \a ptr, данные из \a size элементов типа \a dataType
	void recv(void* ptr, size_t size, MPI_Datatype dataType){
		MPI_Bcast(ptr, size, dataType, senderRank, comm);
	}
};

/**
отправляет подматрицу по заданному средству коммуникации.
\param connector средство коммуникации, которое  следует использовать при обмене.
Представляет переменную типа PeerConnector или BCastConnector.
\param from итератор, указывающий на первую строку подматрицы
\param to итератор, указывающий за последнюю строку подматрицы
\param useBigSends установка в \c true сериализует матрицу перед отправкой для минимизации числа пересылок,
Иначе каждая строка отправляется отдельно.
*/
template <class MPIConnector, class MatrixIterator>
void doSendSubMatrix(MPIConnector connector, MatrixIterator from, MatrixIterator to, bool useBigSends){
	if (useBigSends){
		SerialData data;
		serializeSubMatrix(data,from,to);
		size_t dataSize=data.size();
		connector.send(&dataSize, sizeof(dataSize), MPI_CHAR);
		if (dataSize==0) return;
		connector.send(&data.front(), dataSize, MPI_INT);
	}else{
		int n=to-from;
		connector.send(&n, 1, MPI_INT);
		if (!n) return;
		vector<int> sizes(n);
		int k=0;
		for (MatrixIterator i=from;i!=to;++i,++k){
			sizes[k] = i->size();
		}
		connector.send(&sizes.front(), sizes.size(), MPI_INT);
		for (;from!=to;++from){
			connector.send(&from->front(), INTS_IN_ROW_ELEMENT*from->size(), MPI_INT);
		}
	}
}

/**
Получает строки матрицы по заданному средству коммуникации.
\param connector средство коммуникации, которое  следует использовать при обмене.
Представляет переменную типа PeerConnector или BCastConnector.
\param m матрица, в которую будут дописаны полученные строки.
\param useBigSends установка в \c true означает ожидание получения матрицы в сериализованном виде с последующей десериализаций.
Иначе каждая строка получается отдельно.
*/
template <class MPIConnector, class CMatrix>
int doRecvToMatrix(MPIConnector connector, CMatrix& m, bool useBigSends){
	if (useBigSends){
		size_t dataSize;
		connector.recv(&dataSize, sizeof(dataSize), MPI_CHAR);
		if (dataSize==0) return 0;
		SerialData data(dataSize);
		connector.recv(&data.front(), dataSize, MPI_INT);
		deSerializeToMatrix(data,m);
		return data[0];
	}else{
		int n;
		connector.recv(&n, 1, MPI_INT);
		if (!n) return 0;
		vector<int> sizes(n);
		connector.recv(&sizes.front(), sizes.size(), MPI_INT);
		m.resize(m.size()+n);
		typename CMatrix::iterator i=m.end()-n;
		for (int k=0;k<n;++k,++i){
			i->resize(sizes[k]);
			connector.recv(&i->front(), INTS_IN_ROW_ELEMENT*i->size(), MPI_INT);
		}
		return n;
	}
}


template <class MatrixIterator>
void sendSubMatrix(int recvID, MatrixIterator from, MatrixIterator to, bool useBigSends){
	doSendSubMatrix(PeerConnector(recvID),from,to,useBigSends);
}


template <class CMatrix>
int recvToMatrix(int senderID, CMatrix& m, bool useBigSends){
	return doRecvToMatrix(PeerConnector(senderID),m,useBigSends);
}

template <class MatrixIterator, class MPI_Comm>
void bcastSendSubMatrix(int senderID, MatrixIterator from, MatrixIterator to, MPI_Comm comm, bool useBigSends){
	doSendSubMatrix(BCastConnector(senderID,comm),from,to,useBigSends);
}


template <class CMatrix, class MPI_Comm>
int bcastRecvToMatrix(int senderID, CMatrix& m, MPI_Comm comm, bool useBigSends){
	return doRecvToMatrix(BCastConnector(senderID,comm),m,useBigSends);
}


template <class CMatrix, class MPI_Comm>
void bcastMat(int senderID, CMatrix& m, int myID, MPI_Comm comm, bool useBigSends){
	if (senderID==myID){
		bcastSendSubMatrix(senderID, m.begin(), m.end(), comm, useBigSends);
	}else{
		m.clear();
		bcastRecvToMatrix(senderID, m, comm, useBigSends);
	}
}


//явные инстанциирования
template void sendSubMatrix(int recvID, CMatrix::iterator from, CMatrix::iterator to, bool useBigSends);
template int recvToMatrix(int senderID, CMatrix& m, bool useBigSends);
template void bcastSendSubMatrix(int senderID, CMatrix::iterator from, CMatrix::iterator to, MPI_Comm comm, bool useBigSends);
template int bcastRecvToMatrix(int senderID, CMatrix& m, MPI_Comm comm, bool useBigSends);
template void bcastMat(int senderID, CMatrix& m, int myID, MPI_Comm comm, bool useBigSends);
} //namespace F4MPI
