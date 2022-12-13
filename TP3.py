import statistics
import gzip
import time
import tracemalloc
from queue import PriorityQueue

class LinkedNode:
    def __init__(self,vtx):
        self.vtx = vtx
        self.next = None

class LinkedWeightedNode:
    def __init__(self, vtx, weight=0):
        self.vtx = vtx
        self.next = None
        self.weight = weight

class LinkedFlowNode:
    def __init__(self,vtx,capacity=0,flow=0,reverse_edge = None):
        self.vtx = vtx
        self.next = None
        self.capacity = capacity
        self.flow = flow
        self.reverse_edge = reverse_edge


class ArvoreGeradora:
    def __init__(self,vtx,pai,nivel):
        self.vtx = vtx
        self.pai = pai
        self.nivel = nivel

class Graph:
    def __init__(self, arquivo_relativo ,implementacao=1, directed = False):

        #utilizar apenas with open se for abrir a rede de colaboracao
        with gzip.open(arquivo_relativo, 'rt') as f:
            lines = f.readlines()
            lines = [line.rstrip() for line in lines]
            num_vtx = int(lines[0])
            lines.pop(0)
            for i in range(len(lines)):
                new_edge = lines[i].split()
                new_edge[0] = int(new_edge[0])
                new_edge[1] = int(new_edge[1])
                if(len(new_edge) == 3): #grafo com pesos
                    self.comPeso = True
                    new_edge[2] = float(new_edge[2])
                else:
                    self.comPeso = False
                lines[i] = new_edge

            self.num_vtx = num_vtx
            self.edges = lines
            self.directed = directed

        if(self.comPeso == False):

            self.implementacao = implementacao
            if self.implementacao == 2:
                self.graph_matrix = self.adjacency_matrix()
            if self.implementacao == 1:
                self.graph_list = self.adjacency_list()
            self.graus = self.calcula_graus()
            self.componentes = self.contar_componentes()
            with open('info_grafo.txt', 'w') as f:
                f.write(f'Numero de vertices: {self.num_vtx}// Numero de arestas: {len(self.edges)}// gmin: {self.gmin()}// gmax:{self.gmax()}// N Componentes: {len(self.componentes[0])} // Maior: {max(self.componentes[1])} // Menor: {min(self.componentes[1])}')

        else:

            self.implementacao = 1
            self.graph_list = self.adjacency_list()
            print(f'Vertices: {self.num_vtx}, arestas: {len(self.edges)}')
            
    def adjacency_matrix(self):
        colunas = self.num_vtx
        linhas = self.num_vtx
        matrix = []
        for i in range(linhas):
            linha = []
            for j in range(colunas):
                linha.append(0)
            matrix.append(linha)
        for edge in self.edges:
            linha = edge[0]
            coluna = edge[1]
            matrix[linha-1][coluna-1] = 1 # arestas começam no índice um
            matrix[coluna-1][linha-1] = 1 # matriz é simétrica
        return matrix

    def adjacency_list(self):

        if(self.comPeso == False):
            list = [None] * self.num_vtx
            for i in range(self.num_vtx):
                list[i] = LinkedNode(i+1)
            for edge in self.edges:
                vtx1 = edge[0]
                vtx2 = edge[1]
                
                node = list[vtx1-1]
                while node.next != None:
                    node = node.next
                node.next = LinkedNode(vtx2)

                node = list[vtx2-1]
                while node.next != None:
                    node = node.next
                node.next = LinkedNode(vtx1)
        else:
            if(self.directed == False):
                list = [None] * self.num_vtx
                for i in range(self.num_vtx):
                    list[i] = LinkedWeightedNode(i+1)
                for edge in self.edges:
                    vtx1 = edge[0]
                    vtx2 = edge[1]
                    weight = edge[2]

                    node = list[vtx1-1]
                    while node.next != None:
                        node = node.next
                    node.next = LinkedWeightedNode(vtx2,weight)

                    node = list[vtx2-1]
                    while node.next != None:
                        node = node.next
                    node.next = LinkedWeightedNode(vtx1,weight)
            else:
                list = [None] * self.num_vtx
                for i in range(self.num_vtx):
                    list[i] = LinkedFlowNode(i+1)
                for edge in self.edges:
                    vtx1 = edge[0]
                    vtx2 = edge[1]
                    capacity = edge[2]

                    node = list[vtx1-1]
                    while node.next != None:
                        node = node.next
                    node.next = LinkedFlowNode(vtx2,capacity)

        return list

    def addEdge(self,u,v,capacity,flow=0,reverse_edge=None):
        node = self.graph_list[u-1]
        while node.next != None:
            node = node.next
        node.next = LinkedFlowNode(v,capacity,flow,reverse_edge)
        return node.next

    def removeEdge(self,u,v):
        antecessor = self.graph_list[u-1]
        while antecessor.next.vtx != v:
            antecessor = antecessor.next
        edge = antecessor.next
        if edge.reverse_edge != None:
            edge.reverse_edge.reverse_edge = None
        antecessor.next = edge.next

    def BFS(self,vtx):
        marcado = [False] * self.num_vtx
        fila = []
        arvore_geradora = []
        pais = [None] * self.num_vtx
        niveis = [None] * self.num_vtx

        fila.append(vtx)
        marcado[vtx-1] = True # como os vértices começam no índice 1, precisamos decrementar a posição em marcado
        contador=0

        pais[vtx-1] = 0
        niveis[vtx-1] = 0

        while fila:
            vtx = fila.pop(0)
            contador +=1
            
            #print (vtx, end = " ")
            arvore_geradora.append(ArvoreGeradora(vtx,pais[vtx-1],niveis[vtx-1]))

            #para todo vizinho de vtx w
            #ordena vizinhos
            #se w não marcado,
            #marca e adiciona na fila
            vizinhos_desmarcados = []
            if self.implementacao == 1:
                node = self.graph_list[vtx-1]
                while node.next != None:
                    node = node.next
                    if marcado[node.vtx - 1] == False:
                        vizinhos_desmarcados.append(node.vtx)
            if self.implementacao == 2:
                node = self.graph_matrix[vtx-1]
                for i in range(self.num_vtx):
                    if node[i] == 1:
                        if marcado[i] == False:
                            vizinhos_desmarcados.append(i+1)

            if len(vizinhos_desmarcados) > 0:
                
                vizinhos_desmarcados = sorted(vizinhos_desmarcados)
                for i in vizinhos_desmarcados:
                    fila.append(i)
                    marcado[i-1] = True
                    pais[i-1] = vtx
                    niveis[i-1] = niveis[vtx-1] + 1  
            
        with open('info_BFS.txt', 'w') as b:
            for k in arvore_geradora:
                b.write(f'Vertice: {k.vtx}, pai: {k.pai}, nivel:{k.nivel} //')
        return arvore_geradora, marcado, pais

    def DFS(self,vtx):
        marcado = [False] * self.num_vtx
        pilha = []
        arvore_geradora = []
        pais = [None] * self.num_vtx
        niveis = [None] * self.num_vtx

        pais[vtx-1] = 0
        niveis[vtx-1] = 0
        pilha.append(vtx)

        while pilha:
            vtx = pilha.pop()
            
            if marcado[vtx-1] == False:
                marcado[vtx-1] = True
                
                #print (vtx, end = " ")
                arvore_geradora.append(ArvoreGeradora(vtx,pais[vtx-1],niveis[vtx-1]))

                vizinhos = []
                if self.implementacao == 1:
                    node = self.graph_list[vtx-1]
                    while node.next != None:
                        node = node.next
                        vizinhos.append(node.vtx)
                if self.implementacao == 2:
                    node = self.graph_matrix[vtx-1]
                    for i in range(self.num_vtx):
                        if node[i] == 1:
                            vizinhos.append(i+1)

                if len(vizinhos) > 0:
                    vizinhos = sorted(vizinhos)
                    #loop de tras p frente add pilha
                    for i in reversed(vizinhos):
                        pilha.append(i)
                        pais[i-1] = vtx
                        niveis[i-1] = niveis[vtx-1] + 1
                 

        with open('info_DFS.txt', 'w') as d:
            for k in arvore_geradora:
                d.write(f'Vertice: {k.vtx}, pai: {k.pai}, nivel:{k.nivel} //')
        return arvore_geradora

    def calcula_distancia(self,vtx1,vtx2):
        BFS_vtx1 = self.BFS(vtx1)
        marcado = BFS_vtx1[1]
        if marcado[vtx2-1] == False:
            print('Nao existe caminho')
            return None
        else: 
          for i in BFS_vtx1[0] : #arvore_geradora
            if i.vtx == vtx2:
                print(f'A distancia e de {i.nivel} vertices')
                return i.nivel

    def get_caminho(self,fonte,destino):
        BFS_fonte = self.BFS(fonte)
        marcado = BFS_fonte[1]
        if marcado[destino-1] == False:
            return None
        else:
            pais = BFS_fonte[2]
            caminho = []
            while destino != fonte:
                caminho.append(destino)
                destino = pais[destino-1]
            caminho.append(fonte)
            caminho = list(reversed(caminho))
            return caminho
        
    def contar_componentes(self):
        marcado = [False] * self.num_vtx
        vtx_por_componente = []
        vertices = []
        componentes = 0 
        soma_contadores = 0
        while soma_contadores != self.num_vtx:
            for i in range(len(marcado)):
                if marcado[i] == False:
                    vtx = i+1 #soma pq quando marcamos vtx, marcamos na posição vtx-1
                    break
            BFS_vtx = self.BFS(vtx)
            #iterar dentro do vetor de marcação
            contador = 0
            componente_atual = []
            for i in range(len(BFS_vtx[1])):
                if BFS_vtx[1][i] == True:
                    marcado[i] = True
                    contador += 1
                    soma_contadores += 1
                    componente_atual.append(i+1)
            componentes += 1
            vtx_por_componente.append(contador)
            vertices.append(componente_atual)

        print(f'Temos {componentes} componente(s) neste grafo, maior: {max(vtx_por_componente)}, menor: {min(vtx_por_componente)}')
        with open('info_componentes.txt', 'w') as c:
            for i in range(len(vertices)):
                c.write(f'Componente {i+1} : ')
                for j in vertices[i]:
                    c.write(f'Vertice: {j} //')
        return vertices, vtx_por_componente

    def calcula_diametro(self):
        maximos = [0] * self.num_vtx
        for i in range(1,self.num_vtx + 1):
            BFS_vtx = self.BFS(i)
            maximo_vizinho = BFS_vtx[0].pop() #ultimo elemento do vetor arvore_geradora
            maximo_vizinho = maximo_vizinho.nivel
            maximos[i-1] = maximo_vizinho
        diametro = max(maximos)
        print(f'Diametro: {diametro}')
        return diametro

    def diametro_aproximado(self):
        componentes = self.componentes[0]
        distancia_componente = []
        for i in componentes:
            BFS_primeiro = self.BFS(i[0])[0] #vetor da arvore geradora p primeiro elemento descoberto
            BFS_primeiro = BFS_primeiro.pop()
            BFS_primeiro = BFS_primeiro.nivel #distancia do elemento mais distante
            distancia_componente.append(BFS_primeiro)
        
        return max(distancia_componente)


    def calcula_graus(self):
        if self.implementacao == 2:
            graus = []
            for linha in self.graph_matrix:
                soma = 0
                for i in range(self.num_vtx):
                    soma += linha[i]
                graus.append(soma)
            return graus
        if self.implementacao==1:
            graus = []
            for i in range(self.num_vtx):
                vizinhos = 0
                node = self.graph_list[i]
                while node.next != None:
                    node = node.next
                    vizinhos+=1
                graus.append(vizinhos)
            return graus
    
    def num_edges(self):
        soma_graus = 0
        for grau in self.graus:
            soma_graus += grau
        return int(soma_graus / 2)

    def gmin(self):
        minimo = min(self.graus)
        return minimo

    def gmax(self):
        maximo = max(self.graus)
        return maximo

    def gmed(self):
        soma_graus = 0
        for grau in self.graus:
            soma_graus += grau 
        return soma_graus / self.num_vtx

    def mediana_grau(self):
        return statistics.median(self.graus)

    def dijkstra(self, vtx):
        #testa se é possível
        for e in self.edges:
            if e[2] < float(0):
                print("Pesos negativos !")
                return None
        dist = [1e7] * self.num_vtx
        dist[vtx-1] = 0
        explored = [False] * self.num_vtx
        pai = [None] * self.num_vtx
        
        #explored[v-1] = True se v for explorado
        for contador in range(self.num_vtx):
            dist_min = 1e7
            for i in range(1,self.num_vtx+1): 
                if dist[i-1] < dist_min and explored[i-1] == False:
                    u = i #u vértice de dist min em V-S
                    dist_min = dist[u-1]
            explored[u-1] = True
            v = self.graph_list[u-1].next
            while v != None:
                if (dist[v.vtx - 1] > dist[u-1] + v.weight) and explored[v.vtx - 1] == False: 
                    dist[v.vtx - 1] = dist[u-1] + v.weight
                    pai[v.vtx-1] = u
                v = v.next
        return (dist,pai)

    def Dijkstra_heap(self,vtx):
        for e in self.edges:
            if e[2] < float(0):
                print("Pesos negativos !")
                return None
        dist = [1e7] * self.num_vtx
        dist[vtx-1] = 0
        explored = [False] * self.num_vtx
        pai = [None] * self.num_vtx

        heap = PriorityQueue()
        heap.put((0,vtx))

        for contador in range(self.num_vtx):
            (dist_u,u) = heap.get()
            explored[u-1] = True
            v = self.graph_list[u-1].next
            while v != None:
                if (dist[v.vtx - 1] > dist[u-1] + v.weight) and explored[v.vtx - 1] == False:
                    dist[v.vtx - 1] = dist[u-1] + v.weight
                    pai[v.vtx-1] = u
                    heap.put((dist[v.vtx-1],v.vtx))
                v = v.next
        return (dist,pai)


    def caminho_minimo(self,vtx1,vtx2,comHeap=True):
        if comHeap: 
            dijkstra = self.Dijkstra_heap(vtx1)
        else:
            dijkstra = self.dijkstra(vtx1)
        if dijkstra != None: #não possui pesos negativos
            (dist,pai) = dijkstra
            caminho = []
            print_caminho = 'O caminho minimo e: '
            distancia = dist[vtx2-1]
            if(distancia == 1e7):
                print('Não existe caminho')
                return distancia
            
            while vtx2 != None:
                caminho.append(vtx2)
                vtx2 = pai[vtx2-1]
            for i in reversed(caminho):
                print(i)
                print_caminho += (str(i) + ', ')
            print_caminho += f'distancia: {distancia}'
            print(print_caminho)
            return distancia
        else: 
            return dijkstra

    def Prim(self,vtx):
        custo = [1e7] * self.num_vtx
        explored = [False] * self.num_vtx
        ordem = []
        pai = [None] * self.num_vtx
        custo[vtx-1] = 0
        for contador in range(self.num_vtx):
            custo_min = 1e7
            for i in range(1,self.num_vtx+1): 
                if custo[i-1] < custo_min and explored[i-1] == False:
                    u = i #u vértice de custo min em V-S
                    custo_min = custo[u-1]
            explored[u-1] = True
            ordem.append(u)
            v = self.graph_list[u-1].next
            while v != None:
                if (custo[v.vtx-1] > v.weight) and explored[v.vtx - 1] == False:
                    custo[v.vtx-1] = v.weight
                    pai[v.vtx-1] = u
                v = v.next
        return (custo, ordem, pai)
        
    def print_mst(self,vtx):
        mst = self.Prim(vtx)
        custo_total = 0
        for i in mst[0]:
            if i != None and i != 1e7:
                custo_total += i
        with open('info_mst.txt', 'w') as m:
            for i in range(1,self.num_vtx):
                vertice_explorado = mst[1][i]
                m.write(f'Aresta {mst[2][vertice_explorado-1]} - {vertice_explorado} // ')
            m.write(f'Peso da MST: {custo_total}')

    def get_gargalo(self, aug_path):
        residual_capacities = []
        for i in range(len(aug_path)-1):
            u = aug_path[i]
            v = aug_path[i+1]
            edge1 = self.graph_list[u-1]
            while edge1.vtx != v:
                edge1 = edge1.next
            residual_capacity = edge1.capacity - edge1.flow
            residual_capacities.append(residual_capacity)
            if edge1.reverse_edge == None:
                if edge1.flow < 0: #aresta reversa -> original cheia
                    edge1.reverse_edge = self.addEdge(v,u,residual_capacity,residual_capacity,edge1)
                else:
                    edge1.reverse_edge = self.addEdge(v,u,0,0,edge1)
        gargalo = min(residual_capacities)
        return gargalo
        

    def Ford_Fulkerson(self, s=1, t=2):
        maxflow = 0
        G_f = self
        aug_path = G_f.get_caminho(s,t)
        while aug_path != None:
            gargalo = G_f.get_gargalo(aug_path)
            maxflow += gargalo
            for i in range(len(aug_path)-1):
                u = aug_path[i]
                v = aug_path[i+1]
                edge1 = G_f.graph_list[u-1]
                while edge1.vtx != v:
                    edge1 = edge1.next
                edge1.flow += gargalo
                edge1.reverse_edge.flow -= gargalo
                residual_capacity1 = edge1.capacity - edge1.flow
                residual_capacity2 = edge1.reverse_edge.capacity - edge1.reverse_edge.flow
                #remove aresta se diferença = 0
                if residual_capacity1 == 0:
                    G_f.removeEdge(u,v)
                if residual_capacity2 == 0:
                    G_f.removeEdge(v,u)
            aug_path = G_f.get_caminho(s,t)
        print(f'Maxflow: {maxflow}')
        return maxflow, G_f

    def Flow_allocation(self, s=1, t=2):
        FF = self.Ford_Fulkerson(s,t)
        G_f = FF[1]
        for i in range(self.num_vtx):
            edge = G_f.graph_list[i]
            u = i+1
            while edge.next != None:
                v = edge.vtx
                residual_capacity = edge.capacity - edge.flow
                if edge.flow < 0:
                    if edge.reverse_edge == None:
                        edge.reverse_edge = G_f.addEdge(v,u,residual_capacity,residual_capacity,edge)
                    G_f.removeEdge(u,v)
                edge = edge.next
        with open('alocacao_fluxo.txt', 'w') as f:
            for i in range(self.num_vtx):
                edge = G_f.graph_list[i]
                u = i+1
                while edge.next != None:
                    v = edge.vtx
                    flow = edge.flow
                    capacity = edge.capacity
                    f.write(f'Aresta {u} - {v} , capacidade: {capacity} fluxo: {flow} \n')
                    edge = edge.next
        

#tracemalloc.start()
#G4 = Graph('Teoria dos Grafos\grafos_TP1\grafo_3.txt.gz',1)
#print(tracemalloc.get_traced_memory())
#tracemalloc.stop()

def mil_BFS(G):
    inicio = time.time()
    for i in range(1000):
        if i%100 == 0:
            print('Ainda to aqui')
        i = i % G.num_vtx
        G.BFS(i)   
    fim = time.time()
    tempo_medio = (fim - inicio) / 1000
    print(f'Tempo de exec: {tempo_medio} ')
    return tempo_medio

def mil_DFS(G):
    inicio = time.time()
    for i in range(1000):
        if i%100 == 0:
            print('Ainda to aqui')
        i = i % G.num_vtx
        G.DFS(i)
    fim = time.time()
    tempo_medio = (fim - inicio) / 1000
    print(f'Tempo de exec: {tempo_medio} ')
    return tempo_medio

def cem_dijkstra(G):
    inicio = time.time()
    for i in range(100):
        if i%10 == 0:
            print('Ainda to aqui')
        G.dijkstra(i)
    fim = time.time()
    tempo_medio = (fim - inicio) / 100
    print(f'Tempo de exec: {tempo_medio} ')
    return tempo_medio

def cem_dijkstra_heap(G):
    inicio = time.time()
    for i in range(100):
        if i%10 == 0:
            print('Ainda to aqui')
        G.Dijkstra_heap(i)
    fim = time.time()
    tempo_medio = (fim - inicio) / 100
    print(f'Tempo medio de exec: {tempo_medio} ')
    return tempo_medio

def dez_ff(G):
    inicio = time.time()
    for i in range(10):
        print('Ainda aqui')
        G.Ford_Fulkerson()
    fim = time.time()
    tempo_medio = (fim - inicio) / 10
    print(f'Tempo medio de exec: {tempo_medio} ')
    return tempo_medio

def tempo_FF(G):
    inicio = time.time()
    G.Ford_Fulkerson()
    fim = time.time()
    tempo = (fim - inicio)
    print(f'Tempo de exec: {tempo} ')
    return tempo

RF = Graph('Teoria dos Grafos\grafos_TP1\grafo_rf_4.txt.gz',1,True)
#RF.Flow_allocation()
#tempo_FF(RF)

