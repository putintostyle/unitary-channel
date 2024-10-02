% Define the Hadamard gate
H = (1/sqrt(2)) * [1 1; 1 -1];

% Identity matrix for qubit 2 and qubit 3
I2 = eye(2);

% Tensor product for Hadamard on q1
H1 = kron(kron(H, I2), I2); % Hadamard on q1

% Define the CNOT gate for q1 (control) and q2 (target)
CNOT12 = [1 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0;
          0 0 0 1 0 0 0 0;
          0 0 0 0 1 0 0 0;
          0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 0 1;
          0 0 0 0 0 0 1 0];
      
% Define another Hadamard gate acting on q2
H2 = kron(I2, kron(H, I2)); % Hadamard on q2

% Define the CNOT gate for q2 (control) and q3 (target)
CNOT23 = [1 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0;
          0 0 0 1 0 0 0 0;
          0 0 0 0 1 0 0 0;
          0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 0 1;
          0 0 0 0 0 0 1 0];
      
% Define the final CNOT gate for q1 (control) and q3 (target)
CNOT13 = [1 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0;
          0 0 0 1 0 0 0 0;
          0 0 0 0 1 0 0 0;
          0 0 0 0 0 1 0 0;
          0 0 0 0 0 0 1 0;
          0 0 0 0 0 0 1 0];
      
% Construct the full matrix representation by multiplying all gates
U = CNOT13 * CNOT23 * H2 * CNOT12 * H1;

% Display the matrix
disp(U);