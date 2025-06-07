unit Slover;

interface

uses
  System.SysUtils, System.Math;

const
  MatrixSize = 3;

Type
  TMatrix3x3 = array[0..2, 0..2] of Double;
  TVector3 = array[0..2] of Double;

  TNode = record
    Id: Integer;
    X: Double;
    Y: Double;
  end;

  TElement = record
    Id: Integer;
    Node1, Node2, Node3: TNode;
    Kx, Ky, Thickness: Double;
  end;

  TSolver = class
  private
    procedure GetElementConductionMatrix(const Elem: TElement; const Area: Double;
      out Mat: TMatrix3x3);
    procedure GetElementConvectionMatrix(H, Area, Thickness: Double;
      out Mat: TMatrix3x3);
    procedure GetElementAmbientTempMatrix(H, Area, Thickness, Tinf: Double;
      out Vec: TVector3);
    procedure GetElementHeatSourceMatrix(Q, Area, Thickness: Double;
      out Vec: TVector3);
  public
    class function TriangleArea(const A, B, C: TNode): Double; static;
    class function LineLength(const A, B: TNode): Double; static;
  end;

implementation

{ TSolver }

procedure TSolver.GetElementConductionMatrix(const Elem: TElement; const Area: Double;
  out Mat: TMatrix3x3);
var
  bi, bj, bk, ci, cj, ck: Double;
  B: array[0..1,0..2] of Double;
  D: array[0..1,0..1] of Double;
  i,j: Integer;
  Volume: Double;
begin
  bi := Elem.Node2.Y - Elem.Node3.Y;
  bj := Elem.Node3.Y - Elem.Node1.Y;
  bk := Elem.Node1.Y - Elem.Node2.Y;

  ci := Elem.Node3.X - Elem.Node2.X;
  cj := Elem.Node1.X - Elem.Node3.X;
  ck := Elem.Node2.X - Elem.Node1.X;

  B[0,0] := (1/(2*Area))*bi; B[0,1] := (1/(2*Area))*bj; B[0,2] := (1/(2*Area))*bk;
  B[1,0] := (1/(2*Area))*ci; B[1,1] := (1/(2*Area))*cj; B[1,2] := (1/(2*Area))*ck;

  D[0,0] := Elem.Kx; D[0,1] := 0;
  D[1,0] := 0;       D[1,1] := Elem.Ky;

  Volume := Area * Elem.Thickness;
  for i := 0 to 2 do
    for j := 0 to 2 do
    begin
      Mat[i,j] := Volume*(B[0,i]*D[0,0]*B[0,j] + B[0,i]*D[0,1]*B[1,j] +
                          B[1,i]*D[1,0]*B[0,j] + B[1,i]*D[1,1]*B[1,j]);
    end;
end;

procedure TSolver.GetElementConvectionMatrix(H, Area, Thickness: Double; out Mat: TMatrix3x3);
var
  C: Double;
  i,j: Integer;
begin
  C := (-2*H*Area)/(12*Thickness);
  for i := 0 to 2 do
    for j := 0 to 2 do
      if i=j then
        Mat[i,j] := C*2
      else
        Mat[i,j] := C*1;
end;

procedure TSolver.GetElementAmbientTempMatrix(H, Area, Thickness, Tinf: Double;
  out Vec: TVector3);
var
  C: Double;
  i: Integer;
begin
  C := (-2*H*Tinf*Area)/(3*Thickness);
  for i := 0 to 2 do
    Vec[i] := C;
end;

procedure TSolver.GetElementHeatSourceMatrix(Q, Area, Thickness: Double; out Vec: TVector3);
var
  C: Double;
  i: Integer;
begin
  C := (Q*Area*Thickness)/3;
  for i := 0 to 2 do
    Vec[i] := C;
end;

class function TSolver.LineLength(const A, B: TNode): Double;
begin
  Result := Sqrt(Sqr(B.X - A.X) + Sqr(B.Y - A.Y));
end;

class function TSolver.TriangleArea(const A, B, C: TNode): Double;
begin
  Result := Abs((A.X*(B.Y-C.Y) + B.X*(C.Y-A.Y) + C.X*(A.Y-B.Y))/2);
end;


procedure TSolver.GetEdgeHeatSourceMatrix(Q1, Q2, Q3, L1, L2, L3, Thickness: Double; out Vec: TVector3);
var
  C1, C2, C3: Double;
begin
  C1 := (Q1*L1*Thickness)/2;
  C2 := (Q2*L2*Thickness)/2;
  C3 := (Q3*L3*Thickness)/2;
  Vec[0] := C1 + C3;
  Vec[1] := C1 + C2;
  Vec[2] := C2 + C3;
end;

procedure TSolver.GetEdgeHeatConvectionMatrix(H1,H2,H3,T1,T2,T3,L1,L2,L3,Thickness: Double; out Vec: TVector3);
var
  C1,C2,C3: Double;
begin
  C1 := (H1*T1*L1*Thickness)/2;
  C2 := (H2*T2*L2*Thickness)/2;
  C3 := (H3*T3*L3*Thickness)/2;
  Vec[0] := C1 + C3;
  Vec[1] := C1 + C2;
  Vec[2] := C2 + C3;
end;

procedure TSolver.GetEdgeSpecTempMatrix(T1,T2,T3: Double; out Vec: TVector3);
begin
  Vec[0] := (T1+T3)/2;
  Vec[1] := (T1+T2)/2;
  Vec[2] := (T2+T3)/2;
end;

end.
