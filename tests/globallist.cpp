#include "misc/container.h"

using namespace DROPS;
using namespace std;

typedef NewGlobalListCL<int> iList;



void
Info(iList& l)
{
    cout << "size(): " << l.size()
         << "\tGetNumLevel(): " << l.GetNumLevel()
         << "\tIsEmpty(): " << l.IsEmpty()
         << endl;
}

template<class itT>
void
SeqOut(itT b, itT e)
{
    for (itT it= b; it != e; ++it)
        cout << *it << '\t';
    cout << endl;
}

int main()
{
    typedef iList::iterator itT;
    vector<int*> ad, ad2;
    iList l;
    list<int> l1, l2;
    int i= 0;

    Info( l);
    for (int i= 0; i < 4; ++i) {
        l.AppendLevel();
        Info( l);
        for (int j= 0; j < 3; ++j) {
            l[i].push_back( i*3 + j);
            ad.push_back( &l[i].back());
        }
        cout << "IsLevelEmpty(): " << l.IsLevelEmpty( i) << endl;
        SeqOut( l[i].begin(), l[i].end());
    }
    l.FinalizeModify();

    cout << "Ueberpruefe Adressen auf Gleichheit:" << endl;
    i= 0;
    for (itT it= l.begin(); it != l.end(); ++it, ++i)
        if ( &*it != ad[i] )
            cout << " Adresse verschieden fuer: " << i << endl;
        else
            cout << " Adresse gleich fuer:      " << i << endl;

    Info( l);
    SeqOut( l.begin(), l.end());        
    for (int i= 0; i < 4; ++i)
        SeqOut( l.level_begin( i), l.level_end( i));

    l.PrepareModify();
    Info( l);
    for (int i= 0; i < 4; ++i)
        SeqOut( l[i].begin(), l[i].end());
    l[0].push_back( 22);
    ad.insert( ad.begin() + 3, &l[0].back());
    SeqOut( l[0].begin(), l[0].end());
    l[3].clear();
    SeqOut( l[3].begin(), l[3].end());
    l.FinalizeModify();

    Info( l);
    SeqOut( l.begin(), l.end());        
    for (int i= 0; i < 4; ++i)
        SeqOut( l.level_begin( i), l.level_end( i));
    cout << "IsLevelEmpty( 0): " << l.IsLevelEmpty( 0) << endl;
    cout << "IsLevelEmpty( 3): " << l.IsLevelEmpty( 3) << endl;
    
    l.PrepareModify();
    l.RemoveLastLevel();
    l.FinalizeModify();

    Info( l);
    SeqOut( l.begin(), l.end());        
    for (int i= 0; i < 3; ++i)
        SeqOut( l.level_begin( i), l.level_end( i));

    cout << "Ueberpruefe Adressen auf Gleichheit:" << endl;
    i= 0;
    for (itT it= l.begin(); it != l.end(); ++it, ++i)
        if ( &*it != ad[i] )
            cout << " Adresse verschieden fuer: " << i << endl;

    return 0;
}
