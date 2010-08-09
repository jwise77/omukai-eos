from mpi4py import MPI
import enzo_module as em
gd = em.global_data

class constants:
    shared_state = {}
    def __init__(self):
        self.__dict__ = self.shared_state

c = constants()
for i in dir(em):
    if not i.startswith("E_"): continue
    setattr(c, i[2:], getattr(em, i))

def lgrids(start_grid):
    temp = start_grid
    while temp != None:
        yield start_grid
        temp = start_grid.NextGridThisLevel

def non_block():
    gd.CommunicationDirection = c.COMMUNICATION_POST_RECEIVE
    yield 1
    gd.CommunicationDirection = c.COMMUNICATION_SEND
    yield 2
    em.CommunicationReceiveHandler()

def EvolveHierarchy(top_grid, meta_data, exterior, level_array, initial_dt):
    Stop = False
    if meta_data.Time >= meta_data.StopTime: Stop = True
    if meta_data.CycleNumber >= meta_data.StopCycle: Stop = True
    meta_data.StartCPUTime = meta_data.CPUTime = LastCPUTime = MPI.Wtime()
    meta_data.LastCycleCPUTime = 0.0

    # Optimized CTP would go here

    # Add random forcing if need be
    if gd.RandomForcing:
        for g in lgrids(level_array[0]):
            g.GridData.AppendForcingToBaryonFields()
        Exterior.AppendForcingToBaryonFields()

    # Set top grid boundary condition
    gd.CommunicationReceiveIndex = 0
    gd.CommunicationReceiveCurrentDependsOn = c.COMMUNICATION_NO_DEPENDENCE
    
    for cd in non_block():
        for g in lgrids(level_array[0]):
            i = g.GridData.SetExternalBoundaryValues(exterior)
            if i == c.FAIL: exterior.Prepare(g.GridData)
            em.CopyOverlappingZones(g.GridData, meta_data, level_array, 0)

    # Remove forcing field
    if gd.RandomForcing:
        for g in lgrids(level_array[0]):
            g.GridData.DetachForcingFromBaryonFields()
        Exterior.DetachForcingFromBaryonFields()

    #  XXX
    #  Stuff goes in here!  Important stuff!
    #  But I'll come back to it.
    #  XXX

    initial_dt = 0.1
    print "Evolving"
    em.EvolveLevel(meta_data, level_array, 0, initial_dt, exterior)

def EvolveLevel(meta_data, level_array, level, dt_upper, exterior):
    pass

def main(restart, fn):
    top_grid = em.HierarchyEntry()
    meta_data = em.TopGridData()
    exterior = em.ExternalBoundary()
    level_array = em.LevelHierarchyArray()
    em.CommunicationInitialize([])
    #em.run_enzo_main(["-d", "CollapseTest.enzo"])
    gd.debug = 1
    retval = em.SetDefaultGlobalValues(meta_data)
    initial_dt = 0.0

    if restart:
        em.Group_ReadAllData(fn, top_grid, meta_data, exterior)
        em.CommunicationPartitionGrid(top_grid, 0)
        #gd.CommunicationDirection = c.COMMUNICATION_SEND_RECEIVE
    else:
        initial_dt = em.InitializeNew(fn, top_grid, meta_data, exterior)

    em.AddLevel(level_array, top_grid, 0)
    #em.EvolveHierarchy(top_grid, meta_data, exterior, level_array, initial_dt)
    EvolveHierarchy(top_grid, meta_data, exterior, level_array, initial_dt)

    return top_grid, meta_data, exterior, level_array, initial_dt